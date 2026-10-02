function [theta, Hm0, Tm01, ice_diag, pd] = run_sharc_dsp(T, label, HEAVE_IS_DISPLACEMENT, SCALE_HEAVE, SIGN_ANGLE_DEG, SIGN_QUAD, alpha_c, opts)
% SHARC DSP chain, matched to sharc_meas.m (synthetic), plus heading rotation.
%
%   tilt (comp. filter, alpha_c) -> Az = g - (...) -> detrend -> FIR LPF (filtfilt)
%   -> decimate M=32 -> RC HP (tau=100 s) -> Hanning, Nfft=512 -> accel->disp
%   -> auto/cross spectra -> a1,b1,a2,b2 -> theta1, spread
%
% Heading: slopes are rotated into earth frame PER SAMPLE (100 Hz, before the LPF)
% using psi(t) = psi0 + integral(yaw rate - bias). psi0 (and, with >=2 fixes over
% >=5 min, a bias correction) come from GNSS yaw fixes inside the burst.
% Direction is computed in three frames:
%   body : bow = North, body +y = East (theta relative to bow, clockwise)
%   stbd : earth frame assuming body +y = starboard
%   port : earth frame assuming body +y = port
% The batch script's +/-180 candidates cover the remaining pitch/roll sign cases.
%
% SIGN_QUAD : -1 reproduces the synthetic chain
% alpha_c   : omitted -> 0.999 (synthetic). [] -> adaptive_alpha
% opts      : fs_imu, gnss_yaw (table TimeUTC, YawEstimate_deg[, YawConfidence]),
%             yaw_rate_bias_dps, yaw_sign, gnss_min_conf, gnss_max_gap_s,
%             y_starboard, make_plots, verbose, yaw_deg, station
%             DESIGN (defaults = sharc_meas.m): tau, M, lpf_cutoff_hz, lpf_order,
%             Nfft, f_min, f_max_valid, f_max_band,
%             f_min_mode: 'fixed' (f_min) or 'auto' (spectral minimum between f_min
%             and 0.8*fp -- cuts the 1/f^2 double-integration noise per burst)
%             REGIME: regime_aware, regime_thr (n_ice, n_open, hf_open, kurt_miz),
%             f_max_band_open, despike, dir_estimator ('mem' or 'dfs', display + Pdir_MEM)
%
% Regime detection uses SHARC data only:
%   tail_n   high-frequency slope of Czz ~ f^-n above 1.3*fp (open water ~4-5; ice is
%            much steeper because attenuation grows like f^2-f^4)
%   hf_frac  share of m0 above 1.5*fp (wind sea present -> open water)
%   kurt_hf  kurtosis of heave accel above 2 Hz (impulsive floe collisions -> loose ice)
% Thresholds are calibration constants: set them once from the batch validation
% (ERA5 SIC), then keep them fixed.
% Regime STATE (causal, as it would run onboard): thresholds shift by 'hyst' toward the
% current regime and a change is accepted only after 'persist' consecutive bursts agree.
% Pass pd.regime_state back in as opts.regime_state for the next burst.
% QC: kurtosis > qc_kurt_max or wave-band yaw rms > qc_yaw_rms_max -> 'invalid'
% (handling/deployment, buoy tumbling); the regime state is not updated by such bursts.
% Heading: yaw_bias_mode 'per_burst' removes the burst's own mean yaw rate (no net
% rotation inside a burst; only wave-frequency yaw is rotated out). 'prior' uses
% yaw_rate_bias_dps, which is only safe if the gyro bias is known and stable.
%#ok<*INUSD>  SIGN_ANGLE_DEG unused

    if nargin < 7, alpha_c = 0.999; end
    if nargin < 8 || isempty(opts), opts = struct(); end
    opts = fill_defaults(opts, struct('fs_imu',100, 'gnss_yaw',[], 'yaw_rate_bias_dps',0, ...
        'yaw_sign',1, 'gnss_min_conf',0, 'gnss_max_gap_s',1200, 'y_starboard',true, ...
        'make_plots',true, 'verbose',true, 'yaw_deg',NaN, 'station','SHARC', ...
        'tau',100, 'M',32, 'lpf_cutoff_hz',1.12, 'lpf_order',64, 'Nfft',512, ...
        'f_min',0.03, 'f_max_valid',0.5, 'f_max_band',0.3, 'f_min_mode','fixed', ...
        'regime_aware',true, 'f_max_band_open',0.4, 'despike',true, 'dir_estimator','mem', ...
        'regime_thr', struct('n_ice',7, 'n_open',5, 'hf_open',0.15, 'kurt_miz',6, 'hyst',1, 'persist',2), ...
        'regime_state', struct('regime',"", 'pending',"", 'count',0), ...
        'fp_min',0.045, 'qc_kurt_max',100, 'qc_yaw_rms_max',0.5, ...
        'yaw_bias_mode','per_burst', 'plot_estimators', {{'dfs','mem'}}));
    opts.regime_thr = fill_defaults(opts.regime_thr, struct('n_ice',7, 'n_open',5, 'hf_open',0.15, ...
        'kurt_miz',6, 'hyst',1, 'persist',2));
    opts.regime_state = fill_defaults(opts.regime_state, struct('regime',"", 'pending',"", 'count',0));

    SCALE_GYRO_RATE = 1 / 65.5;
    g = 9.81;
    fs_imu = opts.fs_imu;
    dt = 1 / fs_imu;

    T = sortrows(T, 'DateTimeUTC');
    n_pts = height(T);
    t_s = (0:n_pts-1)' * dt;

    %% Sample-rate check (timestamps vs assumed fs)
    span_s = seconds(T.DateTimeUTC(end) - T.DateTimeUTC(1));
    fs_est = (n_pts - 1) / max(span_s, eps);
    if abs(fs_est/fs_imu - 1) > 0.01
        warning(['[%s] %d samples over %.0f s of timestamps = %.2f Hz, but the chain assumes %.0f Hz ' ...
            '(%.0f s). Gaps or a different ODR: frequencies scale by fs_true/fs_assumed.'], ...
            label, n_pts, span_s, fs_est, fs_imu, t_s(end));
    end

    %% Channels
    heave_accel = SCALE_HEAVE * T.Heave;
    surge_accel = SCALE_HEAVE * T.Surge;
    sway_accel  = SCALE_HEAVE * T.Sway;
    pitch_rate_dps = SCALE_GYRO_RATE * T.Pitch;
    roll_rate_dps  = SCALE_GYRO_RATE * T.Roll;

    %% Floe-collision indicator: kurtosis of heave accel above 2 Hz (Gaussian = 3)
    [bk, ak] = butter(4, 2/(fs_imu/2), 'high');
    hf = filtfilt(bk, ak, heave_accel - mean(heave_accel));
    kurt_hf = mean(hf.^4) / max(mean(hf.^2)^2, eps);
    n_spikes = 0;
    if opts.despike && kurt_hf > opts.regime_thr.kurt_miz
        [heave_accel, n1] = despike(heave_accel);
        [surge_accel, n2] = despike(surge_accel);
        [sway_accel,  n3] = despike(sway_accel);
        n_spikes = n1 + n2 + n3;
        if opts.verbose, fprintf('[%s] kurtosis %.1f: %d impulsive samples replaced\n', label, kurt_hf, n_spikes); end
    end

    %% Accelerometer tilt (identical expression to sharc_meas.m)
    pitch_accel_deg = atan2(surge_accel, sqrt(sway_accel.^2 + heave_accel.^2)) * 180/pi;
    roll_accel_deg  = atan2( sway_accel, sqrt(surge_accel.^2 + heave_accel.^2)) * 180/pi;

    %% Complementary filter (synthetic recursion, zero initial state)
    if ~isempty(alpha_c)
        alpha_vec = alpha_c * ones(n_pts,1);
    else
        [alpha_vec, ad] = adaptive_alpha(pitch_rate_dps, roll_rate_dps, pitch_accel_deg, roll_accel_deg, fs_imu);
        if opts.verbose, fprintf('[%s] adaptive alpha: mean %.4f\n', label, mean(ad.alpha_k)); end
    end
    pitch_angle_deg = zeros(n_pts,1);
    roll_angle_deg  = zeros(n_pts,1);
    for i = 2:n_pts
        a_i = alpha_vec(i);
        pitch_angle_deg(i) = a_i*(pitch_angle_deg(i-1) + pitch_rate_dps(i)*dt) + (1-a_i)*pitch_accel_deg(i);
        roll_angle_deg(i)  = a_i*(roll_angle_deg(i-1)  + roll_rate_dps(i)*dt)  + (1-a_i)*roll_accel_deg(i);
    end
    Sx_fused = deg2rad(pitch_angle_deg);   % along body x (bow)
    Sy_fused = deg2rad(roll_angle_deg);    % along body y

    %% Tilt-compensated vertical accel (synthetic form and sign)
    th = Sx_fused; ph_ = Sy_fused;
    heave_accel = g - (sin(th).*surge_accel + sin(ph_).*cos(th).*sway_accel + cos(ph_).*cos(th).*heave_accel);

    %% Heading
    [psi, heading] = estimate_heading(T, t_s, opts, SCALE_GYRO_RATE);
    has_hd = heading.mode ~= "none";
    if opts.verbose
        fprintf('[%s] heading: %s, %d fixes in burst, bias %+.4f deg/s\n', label, heading.mode, ...
            heading.n_fix, heading.bias_dps);
    end

    %% DSP parameters (defaults match sharc_meas.m; set through opts)
    M = opts.M;
    fs_dec = fs_imu / M;
    if opts.lpf_cutoff_hz >= fs_dec/2
        warning('[%s] LPF cutoff %.2f Hz >= decimated Nyquist %.2f Hz: aliasing.', label, opts.lpf_cutoff_hz, fs_dec/2);
    end
    lpf = fir1(opts.lpf_order, min(0.9, opts.lpf_cutoff_hz / (fs_imu/2)));
    tau = opts.tau;
    k_rc = exp(-1/(tau*fs_dec));
    f_min = opts.f_min;
    chain = @(x) rc_hp(downsample(filtfilt(lpf, 1, x), M), k_rc);

    heave_hp = rc_hp(downsample(filtfilt(lpf, 1, detrend(heave_accel)), M), k_rc);

    Nfft = opts.Nfft;
    if numel(heave_hp) < Nfft
        error('[%s] Only %d decimated samples; need >= %d.', label, numel(heave_hp), Nfft);
    end
    win = hanning(Nfft);
    freq_axis = (0:Nfft/2-1) * (fs_dec/Nfft);
    valid = freq_axis >= f_min & freq_axis <= opts.f_max_valid;
    band  = freq_axis >= f_min & freq_axis <= opts.f_max_band;
    dsp = @(ex, ey) dir_spectra(heave_hp, ex, ey, win, Nfft, fs_dec, valid, freq_axis, ...
                                HEAVE_IS_DISPLACEMENT, SIGN_QUAD, band);

    %% Body frame: bow = North, +y = East
    Sx_hp = chain(Sx_fused);
    Sy_hp = chain(Sy_fused);
    Db = dsp(Sy_hp, Sx_hp);

    %% Earth frames (rotation per sample, before LPF)
    Ds = []; Dp = [];
    if has_hd
        c = cosd(psi); s = sind(psi);
        Ds = dsp(chain(Sx_fused.*s + Sy_fused.*c), chain(Sx_fused.*c - Sy_fused.*s));   % +y starboard
        Dp = dsp(chain(Sx_fused.*s - Sy_fused.*c), chain(Sx_fused.*c + Sy_fused.*s));   % +y port
    end

    %% Low-frequency cut: fixed f_min, or the spectral minimum between f_min and 0.8*fp
    Czz = Db.B.Czz;
    fa  = freq_axis(:);
    fpb = band(:) & fa >= opts.fp_min;                 % no swell peak below fp_min (~22 s)
    [~, kp] = max(Czz .* fpb);
    fp = freq_axis(kp);
    lowfrac = trapz(fa, Czz .* (fa < 0.05)) / trapz(fa, Czz);   % share of m0 below 0.05 Hz (fixed band)
    f_lo = f_min;
    if strcmpi(opts.f_min_mode, 'auto')
        sr = find(fa >= f_min & fa <= 0.8*fp);
        if numel(sr) >= 3
            [~, j] = min(Czz(sr));
            f_lo = fa(sr(j));
        end
    end
    [~, kp] = max(Czz .* (fpb & fa >= f_lo));
    fp = freq_axis(kp);
    band = freq_axis >= f_lo & freq_axis <= opts.f_max_band;
    Db = resummarize(Db, band);
    if has_hd, Ds = resummarize(Ds, band); Dp = resummarize(Dp, band); end

    %% Scalars (moments over [f_lo, f_max_valid])
    mom = fa >= f_lo & fa <= opts.f_max_valid;
    m0 = trapz(fa, Czz .* mom);
    m1 = trapz(fa, fa .* Czz .* mom);
    Hm0 = 4*sqrt(m0);
    Tm01 = m0 / m1;

    %% Gyro yaw rotation over the burst (net, after bias prior) -- heading diagnostic
    gyro_dyaw = NaN;
    if ismember('Yaw', T.Properties.VariableNames)
        yr = opts.yaw_sign * SCALE_GYRO_RATE * double(T.Yaw(:)) - opts.yaw_rate_bias_dps;
        yr(isnan(yr)) = 0;
        gyro_dyaw = sum(yr) * dt;
    end

    %% Regime features (SHARC only)
    floor_ = median(Czz(fa >= 0.35 & fa <= opts.f_max_valid));          % displacement noise floor
    hi = fa > 1.3*fp & fa <= 0.3 & Czz > 10*max(floor_, eps);
    tail_n = NaN;
    if nnz(hi) >= 4
        pp = polyfit(log(fa(hi)), log(Czz(hi)), 1);
        tail_n = -pp(1);
    end
    hf_frac = trapz(fa, Czz .* (fa > 1.5*fp & mom)) / max(m0, eps);
    yaw_rms_wave = NaN;
    if ismember('Yaw', T.Properties.VariableNames)
        [bb, ab] = butter(2, [0.05 0.3]/(fs_imu/2));
        yaw_rms_wave = rms(filtfilt(bb, ab, yr));          % deg/s, buoy yawing at wave frequency
    end
    thr = opts.regime_thr;
    qc_ok = kurt_hf <= opts.qc_kurt_max && ~(yaw_rms_wave > opts.qc_yaw_rms_max);
    st = opts.regime_state;
    regime_raw = classify_regime(tail_n, hf_frac, kurt_hf, thr, "");     % no hysteresis
    if ~qc_ok
        regime = "invalid";                                  % state left untouched
        if opts.verbose, fprintf('[%s] QC FAIL (kurt %.0f, yaw rms %.2f deg/s): handling/tumbling\n', ...
                label, kurt_hf, yaw_rms_wave); end
    else
        cand = classify_regime(tail_n, hf_frac, kurt_hf, thr, st.regime);
        if st.regime == "" || cand == st.regime
            st.regime = cand;  st.pending = "";  st.count = 0;
        elseif cand == st.pending
            st.count = st.count + 1;
            if st.count >= thr.persist, st.regime = cand; st.pending = ""; st.count = 0; end
        else
            st.pending = cand;  st.count = 1;
            if thr.persist <= 1, st.regime = cand; st.pending = ""; st.count = 0; end
        end
        regime = st.regime;
    end

    %% Regime-dependent upper band edge: stop where the signal reaches the noise floor
    f_hi = opts.f_max_band;
    if opts.regime_aware
        fcap = opts.f_max_band;
        if regime == "open", fcap = opts.f_max_band_open; end   % keep the wind sea in open water
        % ("invalid" bursts fall through with the default band)
        k_end = find(fa > fp & fa <= fcap & Czz < max(3*floor_, 1e-3*max(Czz)), 1);
        f_hi = fcap;
        if ~isempty(k_end), f_hi = fa(k_end); end
        band = freq_axis >= f_lo & freq_axis <= f_hi;
        Db = resummarize(Db, band);
        if has_hd, Ds = resummarize(Ds, band); Dp = resummarize(Dp, band); end
    end

    [Pss, fw] = pwelch([Sx_hp(:), Sy_hp(:)], hanning(Nfft), Nfft/2, Nfft, fs_dec);
    Pss = sum(Pss, 2);
    slope_snr = mean(Pss(fw >= 0.05 & fw <= 0.3)) / mean(Pss(fw >= 0.6 & fw <= 1.0));

    if opts.verbose
        fprintf(['[%s] regime %s (tail n=%.1f, hf=%.2f, kurt=%.1f) | Hm0 = %.3f m, Tm01 = %.2f s, ' ...
            'fp = %.3f Hz, band %.3f-%.3f Hz, R1 = %.2f, gyro dyaw = %+.0f deg\n'], label, regime, tail_n, ...
            hf_frac, kurt_hf, Hm0, Tm01, fp, f_lo, f_hi, Db.B.R1_sum, gyro_dyaw);
    end

    %% Primary direction output
    if has_hd
        if opts.y_starboard, Dsel = Ds; frame = 'stbd'; else, Dsel = Dp; frame = 'port'; end
    else
        Dsel = Db; frame = 'body';
    end
    theta = Dsel.B.theta_sum;

    %% Uncertainty (statistical) and system partition, in the output frame
    theta_jack_sd = jack_theta(Dsel.B, band);                      % leave-one-segment-out
    nu_eff = 2 * Db.B.nseg * sum(Czz(mom))^2 / max(sum(Czz(mom).^2), eps);   % Welch-Satterthwaite DOF
    hm0_sd = Hm0 * 0.5 * sqrt(2 / nu_eff);
    part = partition_dirs(Dsel.B, fa, band);
    pdir_mem = mem_peak(Dsel.B, peak_bin(Dsel.B.Czz, band));

    %% Plots
    if opts.make_plots
        rot = 0; 
        if frame == "body" && ~isnan(opts.yaw_deg), rot = opts.yaw_deg; frame = 'earth: body + yaw'; end
        info = sprintf('Regime: %s (raw %s, n=%.1f)  \\sigma_{stat} = %.1f%s', regime, regime_raw, tail_n, ...
            theta_jack_sd, char(176));
        plot_direction(Dsel, freq_axis, band, sprintf('%s [%s]', label, frame), rot, opts.station, ...
            opts.plot_estimators, f_lo, info, part);
    end

    %% Attenuation estimator (synthetic band)
    fa = freq_axis(:);
    lowband = fa > f_min & fa < 0.3;
    ice_diag.attenuation_AB = [NaN; NaN]; ice_diag.Hs_hat = NaN; ice_diag.fp_hat = NaN;
    if sum(lowband) >= 4
        try
            [Hs_hat, fp_hat] = fit_PM(fa(lowband), Czz(lowband));
            E_pred = pm_spectrum(Hs_hat, fp_hat, fa);
            fitband = fa > fp_hat & fa < 0.5;
            logT = log(max(Czz(fitband) ./ E_pred(fitband), eps));
            ice_diag.attenuation_AB = [-fa(fitband).^2, -fa(fitband).^4] \ logT;
            ice_diag.Hs_hat = Hs_hat; ice_diag.fp_hat = fp_hat;
        catch ME
            if opts.verbose, fprintf('[%s] PM fit skipped: %s\n', label, ME.message); end
        end
    end

    %% Diagnostics struct for the batch
    pick = @(D, w) iff(isempty(D), NaN, @() D.(w).theta_sum);
    pd = struct();
    pd.fs_est = fs_est;
    pd.fp = fp;
    pd.freq_axis = freq_axis;
    pd.theta_body       = Db.B.theta_sum;  pd.theta_body_w       = Db.W.theta_sum;
    pd.theta_earth_stbd = pick(Ds,'B');    pd.theta_earth_stbd_w = pick(Ds,'W');
    pd.theta_earth_port = pick(Dp,'B');    pd.theta_earth_port_w = pick(Dp,'W');
    pd.R1_summary   = Db.B.R1_sum;         pd.R1_summary_w   = Db.W.R1_sum;   % rotation-invariant
    pd.slope_snr = slope_snr;
    kB = peak_bin(Db.B.Czz, band);  kW = peak_bin(Db.W.Czz, band);
    pd.pdir_body   = Db.B.theta1(kB);   pd.pdir_body_w   = Db.W.theta1(kW);   % peak direction
    pd.spread_peak = Db.B.spread(kB);   pd.spread_peak_w = Db.W.spread(kW);
    pd.nseg = Db.B.nseg;                pd.nseg_w = Db.W.nseg;
    pd.f_lo_Hz = f_lo;
    pd.Tm01_over_Tp = Tm01 * fp;        % ~0.75-0.85 for a clean wind-sea/swell peak; >1 = low-f energy
    pd.m0_frac_below_0p05 = lowfrac;
    pd.gyro_dyaw_deg = gyro_dyaw;
    pd.pdir_earth_stbd = NaN;  pd.pdir_earth_stbd_w = NaN;
    if has_hd
        pd.pdir_earth_stbd   = Ds.B.theta1(peak_bin(Ds.B.Czz, band));
        pd.pdir_earth_stbd_w = Ds.W.theta1(peak_bin(Ds.W.Czz, band));
    end
    pd.settings = struct('alpha_c', alpha_c, 'fc_comp_Hz', (1-alpha_c)/(2*pi*dt), ...
        'tau_s', tau, 'fc_rc_Hz', 1/(2*pi*tau), 'lpf_cutoff_Hz', opts.lpf_cutoff_hz, ...
        'lpf_order', opts.lpf_order, 'M', M, 'fs_dec_Hz', fs_dec, 'Nfft', Nfft, ...
        'df_Hz', fs_dec/Nfft, 'f_min_Hz', f_min, 'f_min_mode', string(opts.f_min_mode), ...
        'f_max_band_Hz', opts.f_max_band, ...
        'SIGN_QUAD', SIGN_QUAD);
    pd.heading = heading;
    if has_hd
        pd.yaw_mean = circmean_deg(psi);
        z = mean(exp(1i*deg2rad(psi)));
        pd.yaw_std = rad2deg(sqrt(-2*log(min(max(abs(z), eps), 1-eps))));
        pu = rad2deg(unwrap(deg2rad(psi)));
        pd.yaw_drift = pu(end) - pu(1);
    else
        pd.yaw_mean = NaN; pd.yaw_std = NaN; pd.yaw_drift = NaN;
    end
    pd.D = struct('body', Db, 'stbd', Ds, 'port', Dp);
    pd.ice_diag = ice_diag;
    pd.regime_raw = regime_raw;  pd.regime_state = st;  pd.qc_ok = qc_ok;
    pd.regime = regime;  pd.tail_n = tail_n;  pd.hf_frac = hf_frac;  pd.kurt_hf = kurt_hf;
    pd.n_spikes = n_spikes;  pd.yaw_rms_wave = yaw_rms_wave;  pd.f_hi_Hz = f_hi;
    pd.theta_jack_sd = theta_jack_sd;  pd.nu_eff = nu_eff;  pd.hm0_sd = hm0_sd;
    pd.part = part;  pd.part_frame = string(frame);  pd.pdir_mem = pdir_mem;
    pd.settings.regime_aware = opts.regime_aware;
    pd.settings.dir_estimator = string(opts.dir_estimator);
    pd.settings.n_ice = thr.n_ice;  pd.settings.n_open = thr.n_open;
    pd.settings.hf_open = thr.hf_open;  pd.settings.kurt_miz = thr.kurt_miz;
    pd.settings.hyst = thr.hyst;  pd.settings.persist = thr.persist;
    pd.settings.yaw_bias_mode = string(opts.yaw_bias_mode);  pd.settings.fp_min = opts.fp_min;
end

%% ===================== local functions =====================

function D = dir_spectra(z, ex, ey, win, Nfft, fs_dec, valid, freq_axis, HID, SIGN_QUAD, band)
% ex = East(-like) slope, ey = North(-like) slope. Bartlett (B) and 50%-overlap Welch (W).
    N = numel(z);
    K_norm = 2 / (fs_dec * sum(win.^2));
    stB = (0:floor(N/Nfft)-1) * Nfft;
    stW = (0:floor((N-Nfft)/(Nfft/2))) * (Nfft/2);
    D.B = avg_spectra(z, ex, ey, stB, win, Nfft, valid, freq_axis, HID, SIGN_QUAD, K_norm, band);
    D.W = avg_spectra(z, ex, ey, stW, win, Nfft, valid, freq_axis, HID, SIGN_QUAD, K_norm, band);
end

function S = avg_spectra(z, x, y, starts, win, Nfft, valid, freq_axis, HID, SIGN_QUAD, K_norm, band)
    nS = numel(starts);
    Z = zeros(Nfft/2, nS); X = Z; Y = Z;
    for s = 1:nS
        idx = starts(s) + (1:Nfft);
        [Z(:,s), X(:,s), Y(:,s)] = fft_segment(z(idx), x(idx), y(idx), win, Nfft, valid, freq_axis, HID);
    end
    S.nseg = nS;
    S.seg = struct('zz', K_norm*abs(Z).^2, 'xx', K_norm*abs(X).^2, 'yy', K_norm*abs(Y).^2, ...
        'qzx', SIGN_QUAD*K_norm*imag(Z.*conj(X)), 'qzy', SIGN_QUAD*K_norm*imag(Z.*conj(Y)), ...
        'cxy', K_norm*real(X.*conj(Y)));
    S.Czz = K_norm * mean(abs(Z).^2, 2);
    S.Cxx = K_norm * mean(abs(X).^2, 2);
    S.Cyy = K_norm * mean(abs(Y).^2, 2);
    Qzx = SIGN_QUAD * K_norm * mean(imag(Z .* conj(X)), 2);
    Qzy = SIGN_QUAD * K_norm * mean(imag(Z .* conj(Y)), 2);
    Cxy = K_norm * mean(real(X .* conj(Y)), 2);
    den = sqrt(S.Czz .* (S.Cxx + S.Cyy));
    S.a1 = Qzx ./ den;
    S.b1 = Qzy ./ den;
    S.a2 = (S.Cxx - S.Cyy) ./ (S.Cxx + S.Cyy);
    S.b2 = 2*Cxy ./ (S.Cxx + S.Cyy);
    S.theta1 = mod(90 - rad2deg(atan2(S.b1, S.a1)), 360);
    r1 = sqrt(S.a1.^2 + S.b1.^2);
    S.spread = rad2deg(sqrt(max(2*(1 - r1), 0)));
    S = summarize(S, band);
end

function S = summarize(S, band)
% Czz-weighted circular mean direction and mean R1 over the band.
    r1 = sqrt(S.a1.^2 + S.b1.^2);
    w = S.Czz(band); w = w / sum(w);
    S.theta_sum = mod(rad2deg(angle(sum(w .* exp(1i*deg2rad(S.theta1(band))), 'omitnan'))), 360);
    S.R1_sum = sum(w .* r1(band), 'omitnan');
end

function D = resummarize(D, band)
    D.B = summarize(D.B, band);
    D.W = summarize(D.W, band);
end

function [Zd, X, Y] = fft_segment(z, sx, sy, win, Nfft, valid, freq_axis, HID)
% Identical to the per-segment block in sharc_meas.m.
    z = z(:); sx = sx(:); sy = sy(:); win = win(:);
    valid = valid(:); f = freq_axis(:);
    Z = fft(z .* win, Nfft); Z = Z(1:Nfft/2);
    Zd = zeros(Nfft/2, 1);
    if HID
        Zd(valid) = Z(valid);
    else
        Zd(valid) = Z(valid) .* (-1 ./ (2*pi*f(valid)).^2);
    end
    X = fft(sx .* win, Nfft); X = X(1:Nfft/2); X(~valid) = 0;
    Y = fft(sy .* win, Nfft); Y = Y(1:Nfft/2); Y(~valid) = 0;
end

function y = rc_hp(x, k)
    y = zeros(size(x));
    for n = 2:numel(x)
        y(n) = k*y(n-1) + k*(x(n) - x(n-1));
    end
end

function [psi, H] = estimate_heading(T, t_s, opts, scale)
% psi(t) = psi0 + cumtrapz(yaw_rate - b)*dt, anchored to GNSS yaw fixes.
    n = numel(t_s);
    psi = nan(n, 1);
    H = struct('mode', "none", 'n_fix', 0, 'bias_dps', NaN);
    G = opts.gnss_yaw;
    if isempty(G) || ~ismember('Yaw', T.Properties.VariableNames), return; end
    if ismember('YawConfidence', G.Properties.VariableNames)
        G = G(G.YawConfidence >= opts.gnss_min_conf, :);
    end
    if isempty(G), return; end

    rate = opts.yaw_sign * scale * double(T.Yaw(:));
    rate(isnan(rate)) = 0;
    dt = t_s(2) - t_s(1);
    b = opts.yaw_rate_bias_dps;
    per_burst = strcmpi(opts.yaw_bias_mode, 'per_burst');
    if per_burst, b = mean(rate); end            % zero net rotation inside the burst
    rel = @(b) cumtrapz(rate - b) * dt;
    gr = rel(b);

    tf = seconds(G.TimeUTC - T.DateTimeUTC(1));
    yf = G.YawEstimate_deg;
    inb = tf >= 0 & tf <= t_s(end);
    H.n_fix = nnz(inb);

    if H.n_fix >= 1
        [tk, o] = sort(tf(inb)); yk = yf(inb); yk = yk(o);
        H.mode = "gnss_anchor";
        if ~per_burst && H.n_fix >= 2 && (tk(end) - tk(1)) >= 300
            e = rad2deg(unwrap(deg2rad(yk - interp1(t_s, gr, tk))));
            p = polyfit(tk, e, 1);
            delta = -p(1);                    % residual gyro bias, deg/s
            if abs(delta) < 0.5
                b = b + delta; gr = rel(b);
                H.mode = "gnss_fit";
            end
        end
        psi0 = circmean_deg(yk - interp1(t_s, gr, tk));
    else
        dist = max(-tf, tf - t_s(end));      % distance outside the burst, s
        [gap, k] = min(dist);
        if isempty(gap) || gap > opts.gnss_max_gap_s, return; end
        if tf(k) < 0, psi0 = yf(k) - gr(1); else, psi0 = yf(k) - gr(end); end
        H.mode = "gnss_nearest";
    end
    psi = mod(psi0 + gr, 360);
    H.bias_dps = b;
end

function plot_direction(D, freq_axis, band, ttl, rot, station, estimator, f_lo, info, part)
% rot: degrees added to every direction for display (yaw, when theta is body-frame).
% Heatmap: MEM (Lygre & Krogstad 1986) or truncated Fourier (DFS) D(f,theta).
% Pdir = theta1 at the Czz peak, Mdir = Czz-weighted circular mean (estimator independent).
    deg = char(176);
    dirs = @(S) deal(mod(S.theta1(peak_bin(S.Czz, band)) + rot, 360), mod(S.theta_sum + rot, 360));
    [pB, mB] = dirs(D.B);  [pW, mW] = dirs(D.W);

    figure('Name', ['Direction recovery: ' ttl]);
    subplot(2,1,1);
    plot(freq_axis(band), mod(D.B.theta1(band) + rot, 360), 'o-'); hold on;
    plot(freq_axis(band), mod(D.W.theta1(band) + rot, 360), 'x--');
    yline(mB, 'b:', sprintf('Mdir %.1f%s', mB, deg), 'Interpreter','none');
    if ~isnan(part.f_split), xline(part.f_split, 'k--', 'swell | wind sea'); end
    ylabel('\theta_1(f) (deg)'); legend('Bartlett','Welch', 'Location','best'); grid on;
    title(['Recovered direction vs frequency -- ' ttl], 'Interpreter','none');
    subplot(2,1,2);
    plot(freq_axis(band), D.B.spread(band), 'o-'); hold on;
    plot(freq_axis(band), D.W.spread(band), 'x--');
    xlabel('Frequency (Hz)'); ylabel('Spread (deg)'); legend('Bartlett','Welch'); grid on;

    disp_band = freq_axis >= f_lo & freq_axis <= 0.3;          % show the full swell/wind-sea range
    theta_grid = (0:2:360)';
    fmax = 0.3;
    [Theta, F] = meshgrid(theta_grid + rot, freq_axis(disp_band));
    Xc = F .* sind(Theta); Yc = F .* cosd(Theta);
    est = {'B', 'Bartlett', pB, mB; 'W', 'Welch', pW, mW};
    if ischar(estimator) || isstring(estimator), estimator = cellstr(estimator); end
    for q = 1:numel(estimator)
    for r = 1:2
        S = D.(est{r,1});
        E = S.Czz .* dir_dist(S, theta_grid, estimator{q});
        E = E(disp_band,:) / max(E(disp_band,:), [], 'all');
        figure;
        pcolor(Xc, Yc, E); shading interp; colormap(jet); axis equal; hold on;
        add_polar_grid(fmax); xlim([-fmax fmax]*1.2); ylim([-fmax fmax]*1.2);
        title(sprintf('%s %s spectrum -- %s', est{r,2}, upper(estimator{q}), ttl), 'Interpreter','none');
        colorbar; caxis([0 1]);
        txt = sprintf('Pdir(%s) = %.1f%s\nMdir(%s) = %.1f%s', station, est{r,3}, deg, station, est{r,4}, deg);
        if ~isnan(part.mdir_wind)
            txt = sprintf('%s\nswell %.0f%s / wind sea %.0f%s (%.0f%%)', txt, mod(part.mdir_swell + rot, 360), ...
                deg, mod(part.mdir_wind + rot, 360), deg, 100*part.e_wind);
        end
        text(0.02, 0.98, txt, 'Units','normalized', 'VerticalAlignment','top', ...
            'FontSize', 10, 'Interpreter','none', 'BackgroundColor','w');
        text(0.02, 0.02, info, 'Units','normalized', 'VerticalAlignment','bottom', 'FontSize', 9, ...
            'BackgroundColor','w');
    end
    end
end

function g = classify_regime(tail_n, hf_frac, kurt_hf, thr, current)
% Thresholds shift by thr.hyst in favour of the current regime.
    n_ice = thr.n_ice;  n_open = thr.n_open;  hf_open = thr.hf_open;  h = thr.hyst;
    switch current
        case "ice",  n_ice  = n_ice - h;
        case "open", n_open = n_open + h;  hf_open = 0.7*hf_open;
        case "miz",  n_ice  = n_ice + h;   n_open = n_open - h;
    end
    if kurt_hf > thr.kurt_miz
        g = "miz";                                         % impulsive collisions: loose ice
    elseif ~isnan(tail_n) && tail_n >= n_ice && hf_frac < hf_open
        g = "ice";
    elseif (~isnan(tail_n) && tail_n <= n_open) || hf_frac >= hf_open
        g = "open";
    else
        g = "miz";
    end
end

function D = dir_dist(S, theta_grid, estimator)
% Directional distribution per frequency (rows) on compass theta_grid (deg, columns).
    th = deg2rad(theta_grid(:)');
    if strcmpi(estimator, 'mem')
        D = mem_D(S, pi/2 - th);                 % compass -> math angle
    else
        D = max((1/pi)*(0.5 + S.a1.*sin(th) + S.b1.*cos(th) - S.a2.*cos(2*th) + S.b2.*sin(2*th)), 0);
    end
end

function D = mem_D(S, phi)
% Maximum entropy method (Lygre & Krogstad 1986). In the math angle phi our coefficients
% take the standard form D = (1/pi)(1/2 + a1 cos + b1 sin + a2 cos2 + b2 sin2).
    c1 = S.a1 + 1i*S.b1;  c2 = S.a2 + 1i*S.b2;
    r = abs(c1);  c1(r >= 0.999) = c1(r >= 0.999) ./ r(r >= 0.999) * 0.999;
    p1 = (c1 - c2.*conj(c1)) ./ (1 - abs(c1).^2);
    p2 = c2 - c1.*p1;
    num = real(1 - p1.*conj(c1) - p2.*conj(c2));
    den = abs(1 - p1.*exp(-1i*phi) - p2.*exp(-2i*phi)).^2;
    D = max(num, 0) ./ (2*pi*den);
    D(~isfinite(D)) = 0;
end

function th = mem_peak(S, k)
% Compass direction of the MEM maximum at bin k.
    tg = (0:1:359)';
    Sk = struct('a1',S.a1(k),'b1',S.b1(k),'a2',S.a2(k),'b2',S.b2(k));
    [~, j] = max(mem_D(Sk, pi/2 - deg2rad(tg')));
    th = tg(j);
end

function sd = jack_theta(S, band)
% Jackknife standard error (deg) of the Czz-weighted mean direction over segments.
    n = S.nseg;
    if n < 3, sd = NaN; return; end
    th = nan(n,1);
    for j = 1:n
        k = setdiff(1:n, j);
        zz = mean(S.seg.zz(:,k),2);  xx = mean(S.seg.xx(:,k),2);  yy = mean(S.seg.yy(:,k),2);
        qx = mean(S.seg.qzx(:,k),2); qy = mean(S.seg.qzy(:,k),2);
        den = sqrt(zz.*(xx+yy));
        t1 = mod(90 - rad2deg(atan2(qy./den, qx./den)), 360);
        w = zz(band); w = w/sum(w);
        th(j) = rad2deg(angle(sum(w .* exp(1i*deg2rad(t1(band))), 'omitnan')));
    end
    mu = rad2deg(angle(mean(exp(1i*deg2rad(th)))));
    d = mod(th - mu + 180, 360) - 180;
    sd = sqrt((n-1)/n * sum(d.^2));
end

function P = partition_dirs(S, fa, band)
% Split the band at the spectral valley between the two largest peaks (swell = lower f).
    P = struct('n_peaks',0, 'f_swell',NaN, 'f_wind',NaN, 'f_split',NaN, ...
               'mdir_swell',NaN, 'mdir_wind',NaN, 'e_wind',0);
    c = movmean(S.Czz, 3);  c(~band(:)) = 0;
    if ~any(c), return; end
    [pk, loc] = findpeaks(c, 'MinPeakProminence', 0.1*max(c));
    P.n_peaks = numel(pk);
    wmean = @(m) mod(rad2deg(angle(sum(S.Czz(m) .* exp(1i*deg2rad(S.theta1(m))), 'omitnan'))), 360);
    if numel(pk) < 2
        [~, k] = max(c);  P.f_swell = fa(k);  P.mdir_swell = wmean(band(:));
        return;
    end
    [~, o] = sort(pk, 'descend');  two = sort(loc(o(1:2)));
    [~, kv] = min(c(two(1):two(2)));  kv = two(1) + kv - 1;
    P.f_swell = fa(two(1));  P.f_wind = fa(two(2));  P.f_split = fa(kv);
    ms = band(:) & fa < fa(kv);  mw = band(:) & fa >= fa(kv);
    P.mdir_swell = wmean(ms);  P.mdir_wind = wmean(mw);
    P.e_wind = sum(S.Czz(mw)) / sum(S.Czz(band(:)));
end

function [x, nrep] = despike(x)
% Replace impulsive samples (floe collisions): |x - movmedian| > 6 robust sigma.
    med = movmedian(x, 51);
    r = x - med;
    s = 1.4826 * median(abs(r - median(r)));
    bad = abs(r) > 6*s;
    x(bad) = med(bad);
    nrep = nnz(bad);
end

function k = peak_bin(Czz, band)
    [~, k] = max(Czz(:) .* band(:));
end

function add_polar_grid(rmax)
    th = linspace(0, 2*pi, 200);
    for r = 0:0.05:rmax
        plot(r*cos(th), r*sin(th), 'Color',[0.7 0.7 0.7], 'LineStyle',':');
    end
    for ang = 0:45:315
        plot([0 rmax*sind(ang)], [0 rmax*cosd(ang)], 'Color',[0.7 0.7 0.7], 'LineStyle',':');
    end
    text(0, rmax*1.05, 'North', 'HorizontalAlignment','center');
    text(rmax*1.08, 0, 'East', 'HorizontalAlignment','left');
    text(0, -rmax*1.1, 'South', 'HorizontalAlignment','center');
    text(-rmax*1.08, 0, 'West', 'HorizontalAlignment','right');
end

function [Hs_hat, fp_hat] = fit_PM(f_band, Pxx_band)
    f_band = f_band(:); Pxx_band = Pxx_band(:);
    pm_model = @(p, f) (0.3125 * p(1)^2 * (1/p(2))^4 ./ f.^5) .* exp(-1.25 * (1./(p(2)*f)).^4);
    [~, ipk] = max(Pxx_band);
    p0 = [4*sqrt(trapz(f_band, Pxx_band)), 1/f_band(ipk)];
    resid = @(p) pm_model(p, f_band) - Pxx_band;
    p_hat = lsqnonlin(resid, p0, [0 0], [Inf Inf], optimoptions('lsqnonlin','Display','off'));
    Hs_hat = p_hat(1);
    fp_hat = 1/p_hat(2);
end

function E = pm_spectrum(Hs, fp, f)
    E = (0.3125 * Hs^2 * fp^4 ./ f.^5) .* exp(-1.25 * (fp./f).^4);
end

function m = circmean_deg(a)
    a = a(~isnan(a));
    if isempty(a), m = NaN; return; end
    m = mod(atan2d(mean(sind(a)), mean(cosd(a))), 360);
end

function s = fill_defaults(s, d)
    fn = fieldnames(d);
    for k = 1:numel(fn)
        if ~isfield(s, fn{k})
            s.(fn{k}) = d.(fn{k});
        end
    end
end

function v = iff(c, a, bfun)
    if c, v = a; else, v = bfun(); end
end

function [alpha_vec, d] = adaptive_alpha(pitch_rate_dps, roll_rate_dps, pitch_acc_deg, roll_acc_deg, fs)
    band = [0.05 0.25]; fc_min = 0.02; fc_max = 0.80; coh_min = 0.30; alpha_def = 0.99;
    win_s = 600; step_s = 300; lambda = 0.5; nfft = 2048;
    dt = 1/fs; n = numel(pitch_rate_dps);
    [b,a] = butter(2, 0.02/(fs/2), 'high');
    hp = @(x) filtfilt(b, a, x(:));
    gp = hp(cumtrapz(pitch_rate_dps(:))*dt);  ap = hp(pitch_acc_deg);
    gr = hp(cumtrapz(roll_rate_dps(:))*dt);   ar = hp(roll_acc_deg);
    W = min(round(win_s*fs), n);  S = round(step_s*fs);
    starts = 1:S:max(n-W+1, 1);  K = numel(starts);
    g_hat = nan(K,1); coh = nan(K,1); alpha_k = nan(K,1); centers = nan(K,1);
    for k = 1:K
        ix = starts(k):starts(k)+W-1;  centers(k) = mean(ix);
        w_ = hamming(nfft); ov = nfft/2;
        [Cgg1, f] = cpsd(gp(ix), gp(ix), w_, ov, nfft, fs);
        Cgg2 = cpsd(gr(ix), gr(ix), w_, ov, nfft, fs);
        Caa1 = cpsd(ap(ix), ap(ix), w_, ov, nfft, fs);
        Caa2 = cpsd(ar(ix), ar(ix), w_, ov, nfft, fs);
        Cga1 = cpsd(gp(ix), ap(ix), w_, ov, nfft, fs);
        Cga2 = cpsd(gr(ix), ar(ix), w_, ov, nfft, fs);
        Cgg = real(Cgg1+Cgg2);  Caa = real(Caa1+Caa2);  Cga = Cga1+Cga2;
        m = f >= band(1) & f <= band(2);
        g_hat(k) = sum(real(Cga(m))) / sum(Cgg(m));
        gam2 = abs(Cga).^2 ./ (Cgg.*Caa);
        coh(k) = sum(Cgg(m).*gam2(m)) / sum(Cgg(m));
        if coh(k) < coh_min
            target = log(1 - alpha_def);
        else
            wt = min(max(g_hat(k), 0), 1);
            fc = fc_min * (fc_max/fc_min)^wt;
            target = log(dt/(1/(2*pi*fc) + dt));
        end
        if k == 1, st = target; else, st = st + lambda*(target - st); end
        alpha_k(k) = 1 - exp(st);
    end
    if K == 1
        alpha_vec = alpha_k * ones(n,1);
    else
        alpha_vec = interp1(centers, alpha_k, (1:n)', 'linear', 'extrap');
        alpha_vec = min(max(alpha_vec, min(alpha_k)), max(alpha_k));
    end
    d = struct('g_hat',g_hat, 'coh',coh, 'alpha_k',alpha_k, 'centers_s',centers*dt);
end
