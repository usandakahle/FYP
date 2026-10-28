function [theta, Hm0, Tm01, ice_diag, pd] = run_sharc_jack(T, label, HEAVE_IS_DISPLACEMENT, SCALE_HEAVE, SIGN_ANGLE_DEG, SIGN_QUAD, alpha_c, opts)
% SHARC_SIMPLE_V1 (V1.3) + JACKKNIFE 

% Jackknife, K segments used, theta_full from all K:
%   theta_(-i) = direction with segment i left out (same swell region, same masks)
%   d_i  = wrap180(theta_(-i) - theta_full)
%   bias = (K-1) mean(d_i)          theta_jk = theta_full - bias
%   se   = sqrt((K-1)/K sum (d_i - mean d)^2)

% The jackknife removes the first-order (1/K) bias of the estimator. It does not remove a
% bias that every segment shares (e.g. impacts in all of them, or a heading offset).
% With K = 4 it is noisy

% SHARC_SIMPLE_V1 (V1.3) - simple cross-spectral (Longuet-Higgins) directional chain.
%
%   scale -> remove gyro offsets (burst median) -> [despike if impulsive]
%   -> tilt: 2-state Kalman (PI complementary), kP/kI from the Allan variance
%   -> vertical accel g -> detrend (2nd order)
%   -> anti-alias Kaiser FIR (pass f_max_valid, stop fs_dec - f_max_valid, 60 dB), filtfilt
%   -> decimate M -> RC high-pass (tau) on z, Sx, Sy (identical phase on all three)
%   -> Hann, Nfft, Bartlett -> accel -> displacement -> Czz, Cxx, Cyy, Qzx, Qzy, Cxy
%   -> NOISE MODEL of the same chain -> SNR band -> noise corrected a1, b1, a2, b2
%   -> mean direction, peak direction, spread

% IMU NOISE (bench Allan variance):
%   accel white noise n_a = 593 ug/sqrt(Hz), gyro white noise n_g = 0.0352 (deg/s)/sqrt(Hz)

% The noise is passed through the same chain to give its spectrum in each channel:
%   heave displacement : N_zz(f) = (n_a g)^2 |H_hp(f)|^2 / (2 pi f)^4         [m^2/Hz]
%   slope (each axis)  : N_ss(f) = [ n_g^2/(2 pi f)^2 + (n_a kP/(2 pi f))^2 ] |H_hp(f)|^2   [rad^2/Hz]
%   (gyro path integrated above the crossover; accelerometer tilt noise through the
%    complementary low-pass; n in rad). Noise is uncorrelated between channels, so it
%   biases the auto-spectra but not the cross-spectra.

% Band: bins in [f_min, f_max_band] with Czz/N_zz >= snr_min AND (Cxx+Cyy)/(2 N_ss) >= snr_min,
% the contiguous run that contains the peak. This replaces the fixed f_min / 'auto' cut:
% the low edge is where double integration lifts the noise to the signal, per packet.
% Auto-spectra are noise corrected (C - N, floored at 0) before a1-b2, Hm0 and Tm01,
% so low-SNR bins do not pull r1 down or Hm0 up.

% GYRO OFFSETS are removed: pitch and roll rates have zero mean over a
% burst, so the burst median is the offset (turn-on bias + drift). the Kalman bias state
% takes the residual. Yaw offset is handled in the heading (burst mean + ramp test).

% FUSION SIGN CHECK: if the gyro-integrated tilt and the accelerometer tilt are coherent
% (gamma^2 > 0.5) but opposite (gain < -0.2), the gyro is fused with its sign flipped and
% the slope is flipped back afterwards. Wave band slopes stay on the gyro convention
%  only the low frequency fusion stops fighting itself (pd.fuse_flip).

% Inputs: T table with DateTimeUTC, Surge, Sway, Heave, Roll, Pitch
% [, Yaw]; SCALE_HEAVE*raw -> m/s^2; SIGN_ANGLE_DEG and alpha_c are ignored.
% Heading: opts.gnss_yaw (TimeUTC, YawEstimate_deg) or [].

% Outputs: theta (deg, earth if heading, else body; NaN if no SNR band), Hm0 (m,
% noise corrected), Tm01 (s), ice_diag (platform indicators), pd (everything, incl.
% spectra, noise spectra, a1.b2 and a DFS D(f,theta) for plotting).

%#ok<*INUSD>
    if nargin < 8 || isempty(opts), opts = struct(); end
    opts = fill_defaults(opts, struct( ...
        'fs_imu','auto', 'fs_nominal',100, 'gyro_scale_dps',1/65.5, 'acc_offset_mps2',[0 0 0], ...
        'n_acc_ug',593, 'n_gyro_dps',0.0352, 'kf_kP',0.09425, 'kf_kI',0.001773, ...
        'M',32, 'lpf_atten_db',60, 'tau',37, 'Nfft',512, 'detrend_order',2, ...
        'f_min',0.03, 'f_max_valid',0.5, 'f_max_band',0.35, 'fp_min',0.045, ...
        'snr_min',3, 'noise_correct',true, 'despike',true, 'kurt_despike',3.5, ...
        'sys_ratio',0.25, 'sys_trough',0.5, 'sys_halfpower',0.5, 'r1_valid',0.85, ...
        ... impacts and rotation (V1.2)
        'despike_passes',1, 'seg_reject',true, 'seg_kurt_max',3.5, 'seg_rot_max',20, 'seg_min',4, ...
        'kurt_valid',3.5, ...   
        'dir_fmax',0.3, ...    swell-only direction band
        'jackknife',true, 'sigma_valid',NaN, 'hann_bin_factor',2/3, ...  
        'tilt_check_coh',0.5, 'theta_offset_deg',0, ...
        'gnss_yaw',[], 'yaw_sign',1, 'gnss_yaw_sense',1, 'heading_offset_deg',0, ...
        'gnss_min_conf',0, 'gnss_max_gap_s',1200, 'yaw_rate_bias_dps',0, ...
        'yaw_bias_mode','per_burst_ramp', 'ramp_t_min',3, ...
        'verbose',true));
    g = 9.81;

    %% SAMPLES AND RATES
    T = sortrows(T, 'DateTimeUTC');
    n = height(T);
    tt = notz(T.DateTimeUTC);
    span_s = seconds(tt(end) - tt(1));
    fs_ts = (n - 1) / max(span_s, eps);
    if ischar(opts.fs_imu) || isstring(opts.fs_imu)
        fs = opts.fs_nominal;
        if span_s > 60 && abs(fs_ts/opts.fs_nominal - 1) < 0.1, fs = fs_ts; end
    else
        fs = opts.fs_imu;
    end
    dt = 1/fs;
    t_s = (0:n-1)' * dt;
    gaps = seconds(diff(tt));

    %% CHANNELS AND OFFSETS
    ac = SCALE_HEAVE * [double(T.Surge) double(T.Sway) double(T.Heave)] - opts.acc_offset_mps2(:)';
    surge = ac(:,1);  sway = ac(:,2);  heave = ac(:,3);
    p_rate = opts.gyro_scale_dps * double(T.Pitch);
    r_rate = opts.gyro_scale_dps * double(T.Roll);
    gyro_off = [median(p_rate, 'omitnan') median(r_rate, 'omitnan')];
    p_rate = p_rate - gyro_off(1);  r_rate = r_rate - gyro_off(2);

    %% FLOE IMPACTS
    %  Kurtosis of heave accel above 2 Hz: 3 for Gaussian (wave) motion. Despike runs
    % while kurtosis > kurt_despike, up to despike_passes (default 1: in tests, more passes
    % lowered the variance faster than the 4th moment, so kurtosis ROSE). Despiking cannot
    % bring a packet with many impacts back to 3; the segment screen does that.
    % kurt_hf = before, kurt_after = after despike, kurt_used = data actually averaged.
    yr0 = [];
    if ismember('Yaw', T.Properties.VariableNames)
        yr0 = opts.yaw_sign * opts.gyro_scale_dps * double(T.Yaw(:));  yr0(isnan(yr0)) = 0;
        yr0 = yr0 - median(yr0);
    end
    [bk, ak] = butter(4, min(2/(fs/2), 0.9), 'high');
    kurtf = @(x) kurt_of(filtfilt(bk, ak, x - mean(x)));
    kurt_hf = kurtf(heave);  kurt_after = kurt_hf;
    n_spikes = 0;  n_pass = 0;
    while opts.despike && kurt_after > opts.kurt_despike && n_pass < opts.despike_passes
        [surge, n1] = despike(surge);  [sway, n2] = despike(sway);  [heave, n3] = despike(heave);
        [p_rate, n4] = despike(p_rate);  [r_rate, n5] = despike(r_rate);  n6 = 0;
        if ~isempty(yr0), [yr0, n6] = despike(yr0); end
        n_spikes = n_spikes + n1 + n2 + n3 + n4 + n5 + n6;
        n_pass = n_pass + 1;
        kurt_after = kurtf(heave);
        if n1 + n2 + n3 + n4 + n5 + n6 == 0, break; end
    end
    hf = filtfilt(bk, ak, heave - mean(heave));

    %% TILT (Kalman / PI complementary)
    p_acc = atan2d(surge, sqrt(sway.^2 + heave.^2));
    r_acc = atan2d(sway,  sqrt(surge.^2 + heave.^2));
    [gx, cx] = tilt_gain(p_rate, p_acc, fs);
    [gy, cy] = tilt_gain(r_rate, r_acc, fs);
    sx = 1;  sy = 1;
    if cx > opts.tilt_check_coh && gx < -0.2, sx = -1; end
    if cy > opts.tilt_check_coh && gy < -0.2, sy = -1; end
    [p_ang, bp] = pi_tilt(sx*p_rate, p_acc, dt, opts.kf_kP, opts.kf_kI);
    [r_ang, br] = pi_tilt(sy*r_rate, r_acc, dt, opts.kf_kP, opts.kf_kI);
    L60 = min(n, round(60/dt));
    kf_bias = [sx*mean(bp(end-L60+1:end)) sy*mean(br(end-L60+1:end))];

    %% VERT ACCEL AND SLOPES
    th = deg2rad(p_ang);  ph = deg2rad(r_ang);
    az = g - (sin(th).*surge + sin(ph).*cos(th).*sway + cos(ph).*cos(th).*heave);
    Sx = detrend(sx * th, opts.detrend_order);        % wave-band slope on the gyro convention
    Sy = detrend(sy * ph, opts.detrend_order);
    az = detrend(az, opts.detrend_order);

    %% HEADING
    [psi, H] = estimate_heading(T, t_s, opts);
    has_hd = H.mode ~= "none";

    %% ANTI-ALIAS, DECIMATE AND HIGH PASS
    M = opts.M;  fs_dec = fs/M;
    lpf = design_aa(fs, fs_dec, opts.f_max_valid, opts.lpf_atten_db, label);
    k_rc = exp(-1/(opts.tau*fs_dec));
    chain = @(x) rc_hp(downsample(filtfilt(lpf, 1, x), M), k_rc);
    z_d = chain(az);
    Nfft = opts.Nfft;
    if numel(z_d) < Nfft
        error('[%s] %d decimated samples; need >= %d.', label, numel(z_d), Nfft);
    end
    fa = (0:Nfft/2-1)' * (fs_dec/Nfft);  df = fa(2);
    valid = fa >= opts.f_min & fa <= opts.f_max_valid;

    %% NOISE MODEL
    w = 2*pi*fa;  wd = w/fs_dec;
    Hhp2 = abs(k_rc*(1 - exp(-1i*wd)) ./ (1 - k_rc*exp(-1i*wd))).^2;
    na = opts.n_acc_ug*1e-6*9.80665;                   % m/s^2/sqrt(Hz)
    ng = deg2rad(opts.n_gyro_dps);                     % rad/s/sqrt(Hz)
    nat = opts.n_acc_ug*1e-6;                          % accel tilt noise, rad/sqrt(Hz)
    Nzz = na^2 * Hhp2 ./ max(w, eps).^4;
    Nss = (ng^2 ./ max(w, eps).^2 + (nat*opts.kf_kP).^2 ./ max(w, eps).^2) .* Hhp2;
    Nzz(~valid) = 0;  Nss(~valid) = 0;
    if HEAVE_IS_DISPLACEMENT, Nzz(:) = 0; end

    %% SEGMENT SCREEN (Rotations and impacts)
    % Each Bartlett segment (Nfft decimated samples = Nfft*M raw samples) is tested against
    % two known quantities: kurtosis of its HF heave accel (Gaussian = 3, limit seg_kurt_max)
    % and the yaw angle it sweeps (wave-induced yaw is a few degrees, limit seg_rot_max).
    % Segments that fail are left out of the average; if fewer than seg_min pass, the
    % seg_min cleanest are kept and pd.seg_forced is set.
    nS = floor(numel(z_d)/Nfft);
    seg_kurt = nan(nS,1);  seg_rot = zeros(nS,1);
    yint = [];  if ~isempty(yr0), yint = cumsum(yr0)*dt; end
    for q = 1:nS
        ix = ((q-1)*Nfft*M + 1):min(q*Nfft*M, n);
        seg_kurt(q) = kurt_of(hf(ix));
        if ~isempty(yint), seg_rot(q) = max(yint(ix)) - min(yint(ix)); end
    end
    use = true(nS,1);  seg_forced = false;
    if opts.seg_reject
        use = seg_kurt <= opts.seg_kurt_max & seg_rot <= opts.seg_rot_max;
        if nnz(use) < opts.seg_min
            score = seg_kurt/opts.seg_kurt_max + seg_rot/opts.seg_rot_max;
            [~, o] = sort(score);  use(:) = false;  use(o(1:min(opts.seg_min, nS))) = true;
            seg_forced = true;
        end
    end
    iu = cell2mat(arrayfun(@(q) ((q-1)*Nfft*M + 1):min(q*Nfft*M, n), find(use)', 'UniformOutput', false));
    kurt_used = kurt_of(hf(iu));          % kurtosis of the data actually averaged (target ~3)

    %% SPECTRA (body and earth)
    win = hanning(Nfft);
    Sb = xspec(z_d, chain(Sy), chain(Sx), win, Nfft, fs_dec, valid, fa, HEAVE_IS_DISPLACEMENT, SIGN_QUAD, use);
    Se = [];
    if has_hd
        c = cosd(psi);  s = sind(psi);
        Se = xspec(z_d, chain(Sx.*s + Sy.*c), chain(Sx.*c - Sy.*s), win, Nfft, fs_dec, valid, fa, ...
                   HEAVE_IS_DISPLACEMENT, SIGN_QUAD, use);
    end

    %% SNR band 
    Czz = Sb.Czz;
    snr_z = Czz ./ max(Nzz, eps);
    snr_s = (Sb.Cxx + Sb.Cyy) ./ max(2*Nss, eps);
    if HEAVE_IS_DISPLACEMENT, snr_z(:) = Inf; end
    cand = valid & fa >= opts.f_min & fa <= opts.f_max_band & snr_z >= opts.snr_min & snr_s >= opts.snr_min;
    band = false(size(fa));  fp = NaN;  kp = NaN;
    pk = cand & fa >= opts.fp_min;
    if any(pk)
        [~, kp] = max(Czz .* pk);  fp = fa(kp);
        lo = kp;  while lo > 1 && cand(lo-1), lo = lo - 1; end
        hi = kp;  while hi < numel(fa) && cand(hi+1), hi = hi + 1; end
        band(lo:hi) = true;
    end
    f_lo = NaN;  f_hi = NaN;
    if any(band), f_lo = fa(find(band, 1)); f_hi = fa(find(band, 1, 'last')); end

    %% NOISE CORRECTED MOMENTS
    nc = opts.noise_correct;
    % Direction uses only f <= dir_fmax (swell). Above it (wind sea, impact energy that
    % lifts the SNR band to the cap) bins get zero weight in every direction estimate:
    % band mean, swell system, wind-sea system, peak. Hm0/Tm01 still use the full band.
    band_dir  = band & fa <= opts.dir_fmax;
    if ~any(band_dir), band_dir = band; end            % peak itself above dir_fmax: keep band
    kp_dir = kp;
    if ~isnan(kp) && ~band_dir(kp)
        [~, kp_dir] = max(Czz .* band_dir);
    end
    Db = moments(Sb, Nzz, Nss, nc, band_dir, kp_dir, fa);
    De = [];  if has_hd, De = moments(Se, Nzz, Nss, nc, band_dir, kp_dir, fa); end

    Czc = max(Czz - nc*Nzz, 0);

    %% WAVE SYSTEMS (swell / wind sea) 
    % The direction is taken from ONE system, not the whole band: a wide band in a
    % two-peaked sea mixes swell and wind sea, and the mix moves as their energy ratio
    % changes. Systems: local maxima of the 3-bin smoothed, noise-corrected Czz in the
    % band; a second peak counts if >= sys_ratio of the main peak and the trough between
    % them is <= sys_trough of the smaller one. Each system's direction uses its
    % half-power region (Czz >= sys_halfpower * its peak), split at the trough.
    % Primary = the lower-frequency system (swell): it persists for hours; the wind sea
    % follows the wind and is reported separately.
    SY = systems(Czc, band_dir, fa, opts);                       % swell band only
    regD = @(D, m) region_dir(D, Czc, m);
    [th_sw, r1_sw] = regD(Db, SY.m_sw);  [th_ws, r1_ws] = regD(Db, SY.m_ws);
    if has_hd, [th_sw, r1_sw] = regD(De, SY.m_sw);  [th_ws, r1_ws] = regD(De, SY.m_ws); end
    th_sw_body = regD(Db, SY.m_sw);

    % heave-slope coherence, body frame. raw spectra:
    %   gamma^2(f) = (Czx^2 + Qzx^2 + Czy^2 + Qzy^2) / (Czz (Cxx + Cyy))
    % = 1 for a single linear wave train seen by a platform that follows the surface.
    % Lowered by directional spread, sensor noise, and motion that is not wave-driven
    % (impacts, rotation, floe not following). coh_sys: Czz-weighted mean over the swell region.
    coh_f = (Sb.Czx.^2 + Sb.Qzx.^2 + Sb.Czy.^2 + Sb.Qzy.^2) ./ max(Sb.Czz .* (Sb.Cxx + Sb.Cyy), eps);
    coh_f(~valid) = NaN;
    coh_sys = NaN;
    if any(SY.m_sw), coh_sys = sum(coh_f(SY.m_sw) .* Czc(SY.m_sw)) / max(sum(Czc(SY.m_sw)), eps); end
    if ~isfinite(coh_sys) && any(SY.m_sw), coh_sys = mean(coh_f(SY.m_sw), 'omitnan'); end
    dir_valid = any(band) && nnz(SY.m_sw) >= 2 && r1_sw >= opts.r1_valid;

    % r1 alone is not enough: with impacts in the data, r1 >= 0.85 still gave errors of
    % 17-166 deg (SB04, 22-24 Jul). A direction is valid only if the segments averaged are
    % clean (kurt_used <= kurt_valid) and the screen did not have to force segments in.
    invalid_reason = "";
    if ~any(band) || nnz(SY.m_sw) < 2, invalid_reason = "no band";
    elseif r1_sw < opts.r1_valid,       invalid_reason = "r1 below r1_valid";
    elseif seg_forced,                  invalid_reason = "segments forced (impacts/rotation everywhere)";
    elseif kurt_used > opts.kurt_valid, invalid_reason = "impulsive after segment screen";
    end

    % jackknife over the segments used, in the frame that is reported
    if has_hd, Sj = Se; else, Sj = Sb; end
    jk = struct('theta', NaN, 'bias', NaN, 'se', NaN, 'K', Sj.nseg);
    if opts.jackknife && any(SY.m_sw)
        jk = jack_dir(Sj, Nzz, Nss, nc, SY.m_sw, fa, th_sw);
    end

    %  analytic random error of theta1 (Chapter 3, eq:theta_error). For linear waves the
    % heave-slope coherence equals r1^2; K = segments x independent bins in the swell region
    % (Hann: neighbouring bins correlated, hann_bin_factor ~ 2/3 independent per bin).
    K_eff = Sj.nseg * max(1, opts.hann_bin_factor * nnz(SY.m_sw));
    g2 = min(max(r1_sw, 0)^2, 0.999);
    sigma_theta = rad2deg(sqrt((1 - g2) / (2 * K_eff * max(g2, eps))));
    % optional validity by the size of the random error
    if isempty(char(invalid_reason)) && isfinite(opts.sigma_valid) && sigma_theta > opts.sigma_valid
        invalid_reason = "sigma_theta above sigma_valid";
    end
    dir_valid = isempty(char(invalid_reason));

    %% INDICATORS (for flagging)
    %   kurtosis      : 3 for Gaussian wave motion          -> impulsive if kurt_hf > 3.5
    %   spikes        : P(|x| > 6 sigma) = 2e-9 for Gaussian -> any n_spikes > 10 is not wave motion
    %   tail slope    : wave spectrum ~ f^-4 .. f^-5 above the peak -> > -3 = flattened (broadband)
    %   tail check    : (Cxx+Cyy)/(k^2 Czz) = 1 for linear waves, upper half of band -> outside 0.5-2
    %   band at cap   : f_hi pinned at f_max_band -> energy above the IMU noise right up to the cap
    %   rotation      : yaw swept inside a segment -> > seg_rot_max deg
    tail_n = NaN;  chk_tail = NaN;
    if ~isnan(kp)
        tb = band & fa >= 1.3*fp & snr_z >= 10 & Czc > 0;
        if nnz(tb) >= 4, pq = polyfit(log(fa(tb)), log(Czc(tb)), 1);  tail_n = pq(1); end
        kw = (2*pi*fa).^2 / g;
        Rk = (max(Sb.Cxx - nc*Nss, 0) + max(Sb.Cyy - nc*Nss, 0)) ./ (kw.^2 .* max(Czc, eps));
        tc = band & fa >= 1.3*fp;
        if any(tc), chk_tail = median(Rk(tc)); end
    end
    band_cap = ~isnan(f_hi) && f_hi >= opts.f_max_band - df;

    flo_m = opts.f_min;  if ~isnan(f_lo), flo_m = max(opts.f_min, f_lo); end
    mom = valid & fa >= flo_m;
    m0 = sum(Czc(mom))*df;  m1 = sum(fa(mom).*Czc(mom))*df;
    Hm0 = 4*sqrt(m0);  Tm01 = m0/max(m1, eps);
    Hm0_raw = 4*sqrt(sum(Czz(valid))*df);
    Hm0_noise = 4*sqrt(sum(Nzz(valid))*df);

    %% OUTPUTS
    off = opts.theta_offset_deg;
    rot = @(a) mod(a + off, 360);
    if has_hd, Dsel = De; frame = "earth"; else, Dsel = Db; frame = "body"; end
    theta = rot(th_sw);                  % primary: swell system (see wave systems)
    if ~any(band), theta = NaN; end

    yaw_max = NaN;  gyro_dyaw = NaN;
    if ismember('Yaw', T.Properties.VariableNames)
        yr = opts.yaw_sign * opts.gyro_scale_dps * double(T.Yaw(:)); 
        yr(isnan(yr)) = 0;
        yr = yr - median(yr);
        gyro_dyaw = sum(yr)*dt;
        yaw_max = max(abs(movmedian(yr, 5)));
    end
    psi_v = [NaN NaN NaN];
    if has_hd, psi_v = [psi(1) psi(end) circmean_deg(psi)]; end

    ice_diag = struct('kurt_hf', kurt_hf, 'n_spikes', n_spikes, 'tilt_gain', mean([gx gy], 'omitnan'), ...
        'tilt_coh', mean([cx cy], 'omitnan'), 'yaw_max_dps', yaw_max);

    tg = (0:5:355)';
    pd = struct();
    pd.version = "SHARC_SIMPLE_V1_JACK";
    pd.label = string(label);
    pd.t_start = tt(1);  pd.t_end = tt(end);  pd.t_mid = tt(1) + (tt(end) - tt(1))/2;

    pd.n_samples = n;  pd.duration_s = span_s;  pd.fs = fs;  pd.fs_ts = fs_ts;

    pd.max_gap_s = max([gaps; 0]);  pd.M = M;  pd.fs_dec = fs_dec;  pd.Nfft = Nfft;  pd.nseg = Sb.nseg;

    pd.lpf_taps = numel(lpf);

    pd.gyro_offset_dps = gyro_off;  pd.kf_bias_dps = kf_bias;

    pd.fuse_flip = [sx sy] < 0;  pd.tilt_gain = [gx gy];  pd.tilt_coh = [cx cy];

    pd.kurt_hf = kurt_hf;  pd.n_spikes = n_spikes;

    pd.kurt_after = kurt_after;  pd.despike_passes = n_pass;

    pd.seg_kurt = seg_kurt;  pd.seg_rot = seg_rot;  pd.seg_used = use;   

    pd.nseg_total = nS;  pd.nseg_used = nnz(use);  pd.seg_forced = seg_forced;

    pd.seg_rot_max = max(seg_rot);
    pd.kurt_used = kurt_used;

    pd.tail_n = tail_n;  pd.chk_tail = chk_tail;  pd.band_cap = band_cap; 

    pd.yaw_max = yaw_max;  pd.gyro_dyaw = gyro_dyaw;

    pd.heading_mode = H.mode;  pd.heading_n_fix = H.n_fix;  pd.heading_bias_dps = H.bias_dps;

    pd.psi_start = psi_v(1);  pd.psi_end = psi_v(2);  pd.psi_mid = psi_v(3);

    pd.freq = fa;  pd.Czz = Czz;  pd.Cxx = Sb.Cxx;  pd.Cyy = Sb.Cyy;  pd.Nzz = Nzz;  pd.Nss = Nss;

    pd.snr_z = snr_z;  pd.snr_s = snr_s;  pd.band = band;

    pd.f_lo = f_lo;  pd.f_hi = f_hi;  pd.n_band = nnz(band);  pd.fp = fp;  pd.Tp = 1/fp;

    pd.Hm0 = Hm0;  pd.Hm0_raw = Hm0_raw;  pd.Hm0_noise = Hm0_noise;  pd.Tm01 = Tm01;

    pd.snr_z_peak = NaN;  pd.snr_s_peak = NaN;

    if ~isnan(kp), pd.snr_z_peak = snr_z(kp);  pd.snr_s_peak = snr_s(kp); end

    pd.frame = frame;  pd.theta_offset_deg = off;

    pd.theta = theta;  pd.dir_valid = dir_valid;  pd.r1_valid = opts.r1_valid;
    pd.invalid_reason = invalid_reason;  

    pd.theta_jk = rot(jk.theta);  pd.jk_bias = jk.bias;  pd.jk_se = jk.se;  pd.jk_K = jk.K;  

    pd.sigma_theta = sigma_theta;  pd.K_eff = K_eff;  pd.n_bins_sys = nnz(SY.m_sw);

    pd.r1_sys = r1_sw;  pd.spread_sys = rad2deg(sqrt(max(2*(1 - r1_sw), 0)));

    pd.theta_body_sys = rot(th_sw_body);

    pd.bimodal = SY.bimodal;  pd.f_sw = SY.f_sw;  pd.f_ws = SY.f_ws;  pd.f_trough = SY.f_trough;

    pd.theta_ws = rot(th_ws);  pd.r1_ws = r1_ws;  pd.sys_sw = SY.m_sw;  pd.sys_ws = SY.m_ws;

    pd.Hs_sw = 4*sqrt(sum(Czc(SY.side_sw))*df);  pd.Hs_ws = 4*sqrt(sum(Czc(SY.side_ws))*df);
    
    pd.theta_band = rot(Dsel.theta_mean);       % band-mean Mdir (f <= dir_fmax), for comparison
    pd.dir_fmax = opts.dir_fmax;  pd.n_band_dir = nnz(band_dir);   

    pd.theta_body = rot(Db.theta_mean);  pd.theta_body_peak = rot(Db.theta_peak);
    pd.theta_earth = NaN;  pd.theta_earth_peak = NaN;

    if has_hd, pd.theta_earth = rot(De.theta_mean);  pd.theta_earth_peak = rot(De.theta_peak); end

    pd.theta_peak = rot(Dsel.theta_peak);
    pd.spread_mean = Dsel.spread_mean;  pd.spread_peak = Dsel.spread_peak;

    pd.r1_mean = Dsel.r1_mean;  pd.r1_peak = Dsel.r1_peak;  pd.coh_peak = Dsel.coh_peak;
    pd.check_peak = Dsel.check_peak;

    pd.coh_f = coh_f;  pd.coh_sys = coh_sys;                                        
    pd.a1 = Dsel.a1;  pd.b1 = Dsel.b1;  pd.a2 = Dsel.a2;  pd.b2 = Dsel.b2;
    pd.theta1 = rot(Dsel.theta1);  pd.r1 = Dsel.r1;
    pd.theta_grid = tg;  pd.D = dfs(Dsel, tg - off);   % D(f, theta) on the output grid
    pd.qc_note = "";
    if ~any(band), pd.qc_note = "no bins above the noise (snr_min)"; end

    if opts.verbose
        fprintf(['[%s] fs %.2f Hz | Hm0 %.2f m (raw %.2f, noise %.3f) Tm01 %.1f s fp %.3f Hz | band %.3f-%.3f Hz (%d bins) | ' ...
            '%s theta %.0f (peak %.0f), r1 %.2f, spread %.0f | kurt %.1f, tilt gain %.2f/%.2f%s\n'], ...
            label, fs, Hm0, Hm0_raw, Hm0_noise, Tm01, fp, f_lo, f_hi, nnz(band), frame, theta, pd.theta_peak, ...
            pd.r1_mean, pd.spread_mean, kurt_hf, gx, gy, iff(any(pd.fuse_flip), ' (fusion sign flipped)', @() ''));
    end
end

%% LOCAL FUNCTIONS

function S = xspec(z, ex, ey, win, Nfft, fs_dec, valid, fa, HID, SIGN_QUAD, use)
% Bartlett (non-overlapping) auto and cross spectra. ex = east-like, ey = north-like slope.
% use: logical per segment ; only those segments are averaged.
    N = numel(z);  nS = floor(N/Nfft);
    if nargin < 11 || isempty(use), use = true(nS,1); end
    K = 2 / (fs_dec * sum(win.^2));
    Z = zeros(Nfft/2, nS);  X = Z;  Y = Z;
    f = fa(:);
    for s = 1:nS
        idx = (s-1)*Nfft + (1:Nfft);
        Zs = fft(z(idx).*win, Nfft);  Zs = Zs(1:Nfft/2);
        if ~HID, Zs = Zs .* (-1 ./ max(2*pi*f, eps).^2); end
        Zs(~valid) = 0;
        Xs = fft(ex(idx).*win, Nfft);  Xs = Xs(1:Nfft/2);  Xs(~valid) = 0;
        Ys = fft(ey(idx).*win, Nfft);  Ys = Ys(1:Nfft/2);  Ys(~valid) = 0;
        Z(:,s) = Zs;  X(:,s) = Xs;  Y(:,s) = Ys;
    end
    Z = Z(:, use);  X = X(:, use);  Y = Y(:, use);
    S.nseg = nnz(use);
    S.Zs = Z;  S.Xs = X;  S.Ys = Y;  S.K = K;  S.sq = SIGN_QUAD;     %  per-segment spectra
    S.Czz = K*mean(abs(Z).^2, 2);  S.Cxx = K*mean(abs(X).^2, 2);  S.Cyy = K*mean(abs(Y).^2, 2);
    S.Qzx = SIGN_QUAD*K*mean(imag(Z.*conj(X)), 2);
    S.Qzy = SIGN_QUAD*K*mean(imag(Z.*conj(Y)), 2);
    S.Cxy = K*mean(real(X.*conj(Y)), 2);
    S.Czx = K*mean(real(Z.*conj(X)), 2);  S.Czy = K*mean(real(Z.*conj(Y)), 2);   % co-spectra
end

function D = moments(S, Nzz, Nss, nc, band, kp, fa)
% Longuet-Higgins coefficients with noise corrected auto-spectra.
    Czz = max(S.Czz - nc*Nzz, 0);
    Cxx = max(S.Cxx - nc*Nss, 0);  Cyy = max(S.Cyy - nc*Nss, 0);
    den = sqrt(Czz .* (Cxx + Cyy));  ds = Cxx + Cyy;
    D.a1 = S.Qzx ./ den;  D.b1 = S.Qzy ./ den;
    D.a2 = (Cxx - Cyy) ./ ds;  D.b2 = 2*S.Cxy ./ ds;
    D.theta1 = mod(90 - atan2d(D.b1, D.a1), 360);
    D.r1 = min(hypot(D.a1, D.b1), 1);
    kw = (2*pi*fa).^2 / 9.81;
    chk = ds ./ (kw.^2 .* Czz);
    cohr = (S.Qzx.^2 + S.Qzy.^2) ./ (S.Czz .* (S.Cxx + S.Cyy));     % raw, includes noise
    D.theta_mean = NaN;  D.r1_mean = NaN;  D.spread_mean = NaN;
    D.theta_peak = NaN;  D.r1_peak = NaN;  D.spread_peak = NaN;  D.coh_peak = NaN;  D.check_peak = NaN;
    if any(band)
        % band mean direction from energy weighted first moments (standard buoy Mdir)
        A = sum(D.a1(band).*Czz(band), 'omitnan');  B = sum(D.b1(band).*Czz(band), 'omitnan');
        D.theta_mean = mod(90 - atan2d(B, A), 360);
        D.r1_mean = min(hypot(A, B) / sum(Czz(band)), 1);
        D.spread_mean = rad2deg(sqrt(2*(1 - D.r1_mean)));
    end
    if ~isnan(kp)
        k = max(kp-1, 1):min(kp+1, numel(fa));  k = k(band(k));      % peak +/- 1 bin, in band
        A = sum(D.a1(k).*Czz(k), 'omitnan');  B = sum(D.b1(k).*Czz(k), 'omitnan');
        D.theta_peak = mod(90 - atan2d(B, A), 360);
        D.r1_peak = min(hypot(A, B) / sum(Czz(k)), 1);
        D.spread_peak = rad2deg(sqrt(2*(1 - D.r1_peak)));
        D.coh_peak = cohr(kp);  D.check_peak = chk(kp);
    end
end

function SY = systems(C, band, fa, o)
% Swell / wind-sea partition of the band (see 'wave systems' in the main body).
    n = numel(C);  SY = struct('bimodal', false, 'f_sw', NaN, 'f_ws', NaN, 'f_trough', NaN, ...
        'm_sw', false(n,1), 'm_ws', false(n,1), 'side_sw', false(n,1), 'side_ws', false(n,1));
    kb = find(band);
    if numel(kb) < 2, SY.m_sw = band; SY.side_sw = band; return; end
    cs = C;  cs(kb) = conv(C(kb), ones(3,1)/3, 'same');
    cs(kb([1 end])) = C(kb([1 end]));
    pk = [];
    for j = 1:numel(kb)
        k = kb(j);
        l = -Inf;  if j > 1, l = cs(kb(j-1)); end
        r = -Inf;  if j < numel(kb), r = cs(kb(j+1)); end
        if cs(k) > l && cs(k) >= r, pk(end+1) = k; end %#ok<AGROW>
    end
    [~, im] = max(cs(pk));  km = pk(im);
    k2 = NaN;  best = 0;
    for k = pk
        if k == km || cs(k) < o.sys_ratio*cs(km), continue; end
        seg = min(k, km):max(k, km);
        if min(cs(seg)) <= o.sys_trough*min(cs(k), cs(km)) && cs(k) > best, k2 = k; best = cs(k); end
    end
    side = band;  kt = NaN;
    if ~isnan(k2)
        SY.bimodal = true;
        seg = min(k2, km):max(k2, km);  [~, j] = min(cs(seg));  kt = seg(j);
        SY.f_trough = fa(kt);
        k1 = min(k2, km);  k2 = max(k2, km);
    else
        k1 = km;
    end
    hp = @(kc, sd) grow(cs >= o.sys_halfpower*cs(kc) & sd, kc);
    if SY.bimodal
        sd1 = band & (1:n)' < kt;  sd2 = band & (1:n)' > kt;
        SY.m_sw = hp(k1, sd1);  SY.m_ws = hp(k2, sd2);  SY.f_ws = fa(k2);
        SY.side_sw = sd1;  SY.side_ws = sd2;
    else
        SY.m_sw = hp(k1, side);  SY.side_sw = band;
    end
    SY.f_sw = fa(k1);
end

function m = grow(ok, k)
% Contiguous run of ok containing k.
    m = false(size(ok));
    if ~ok(k), m(k) = true; return; end
    lo = k;  while lo > 1 && ok(lo-1), lo = lo - 1; end
    hi = k;  while hi < numel(ok) && ok(hi+1), hi = hi + 1; end
    m(lo:hi) = true;
end

function [th, r1] = region_dir(D, C, m)
% Energy weighted first-moment direction and r1 over the bins in m.
    th = NaN;  r1 = NaN;
    if ~any(m), return; end
    A = sum(D.a1(m).*C(m), 'omitnan');  B = sum(D.b1(m).*C(m), 'omitnan');
    th = mod(90 - atan2d(B, A), 360);
    r1 = min(hypot(A, B) / max(sum(C(m)), eps), 1);
end

function jk = jack_dir(S, Nzz, Nss, nc, m, fa, th_full)
% Leave one segment out swell direction (same region m, same noise correction).
    K = size(S.Zs, 2);
    jk = struct('theta', NaN, 'bias', NaN, 'se', NaN, 'K', K);
    if K < 3 || isnan(th_full), return; end
    d = nan(K, 1);
    for i = 1:K
        k = true(1, K);  k(i) = false;
        Z = S.Zs(:, k);  X = S.Xs(:, k);  Y = S.Ys(:, k);
        Si.Czz = S.K*mean(abs(Z).^2, 2);  Si.Cxx = S.K*mean(abs(X).^2, 2);  Si.Cyy = S.K*mean(abs(Y).^2, 2);
        Si.Qzx = S.sq*S.K*mean(imag(Z.*conj(X)), 2);  Si.Qzy = S.sq*S.K*mean(imag(Z.*conj(Y)), 2);
        Si.Cxy = S.K*mean(real(X.*conj(Y)), 2);
        Di = moments(Si, Nzz, Nss, nc, m, NaN, fa);
        thi = region_dir(Di, max(Si.Czz - nc*Nzz, 0), m);
        d(i) = mod(thi - th_full + 180, 360) - 180;
    end
    d = d(isfinite(d));  K = numel(d);
    if K < 3, return; end
    jk.bias  = (K - 1) * mean(d);
    jk.theta = mod(th_full - jk.bias, 360);
    jk.se    = sqrt((K - 1)/K * sum((d - mean(d)).^2));
    jk.K     = K;
end

function Dm = dfs(D, tg)
% Weighted truncated Fourier series D(f, theta) (always >= 0 for r1, r2 <= 1).
% Compass angles: theta1 = 90 - atan2(b1, a1), theta2 = 90 - atan2(b2, a2)/2.
    r1 = D.r1(:);  r2 = min(hypot(D.a2(:), D.b2(:)), 1);
    t1 = D.theta1(:);  t2 = 90 - atan2d(D.b2(:), D.a2(:))/2;
    tg = tg(:)';
    Dm = (0.5 + (2/3)*r1.*cosd(tg - t1) + (1/6)*r2.*cosd(2*(tg - t2))) / pi;
    Dm(~isfinite(Dm)) = 0;
end

function h = design_aa(fs, fs_dec, fpass, atten_db, label)
% Kaiser anti-alias FIR: everything that folds onto [0, fpass] after decimation to fs_dec
% (bands k*fs_dec +/- fpass) is attenuated by atten_db, so the white accelerometer and
% gyro noise above fs_dec/2 is not folded into the wave band.
    fstop = fs_dec - fpass;
    if fstop <= fpass
        error('[%s] fs_dec = %.2f Hz too low for f_max_valid = %.2f Hz.', label, fs_dec, fpass);
    end
    [nk, Wn, beta, ftype] = kaiserord([fpass fstop], [1 0], [0.01 10^(-atten_db/20)], fs);
    nk = nk + mod(nk, 2);
    h = fir1(nk, Wn, ftype, kaiser(nk+1, beta), 'noscale');
end

function [th, b] = pi_tilt(rate, acc, dt, kP, kI)
% Steady-state 2-state Kalman (angle, gyro bias) = PI complementary filter, degrees.
    n = numel(rate);  th = zeros(n,1);  b = zeros(n,1);  th(1) = acc(1);
    for i = 2:n
        e = acc(i) - th(i-1);
        b(i)  = b(i-1) - kI*e*dt;
        th(i) = th(i-1) + (rate(i) - b(i) + kP*e)*dt;
    end
end

function y = rc_hp(x, k)
    y = zeros(size(x));
    for i = 2:numel(x), y(i) = k*y(i-1) + k*(x(i) - x(i-1)); end
end

function [gain, coh] = tilt_gain(rate_dps, acc_deg, fs)
% Gain and coherence of accelerometer tilt on gyro-integrated tilt, 0.03-0.15 Hz.
    M = 32;  L = 512;  n = floor(numel(rate_dps)/M);  gain = NaN;  coh = NaN;
    if n < L, return; end
    r = mean(reshape(rate_dps(1:n*M), M, n), 1)';
    a = mean(reshape(acc_deg(1:n*M),  M, n), 1)';
    fsd = fs/M;  w = hanning(L);  hop = L/2;
    f = (0:L-1)'*fsd/L;  bnd = f >= 0.03 & f <= 0.3;  wv = 2*pi*f(bnd);
    Gtt = 0;  Gaa = 0;  Gta = 0;
    for s0 = 1:hop:(n - L + 1)
        rs = r(s0:s0+L-1);  as = a(s0:s0+L-1);
        R = fft((rs - mean(rs)).*w);  A = fft((as - mean(as)).*w);
        Th = R(bnd) ./ (1i*wv);  Ab = A(bnd);
        Gtt = Gtt + sum(abs(Th).^2);  Gaa = Gaa + sum(abs(Ab).^2);  Gta = Gta + sum(conj(Th).*Ab);
    end
    gain = real(Gta)/Gtt;  coh = abs(Gta)^2/(Gtt*Gaa);
end

function k = kurt_of(x)
% kurtosis (3 = Gaussian).
    x = x - mean(x);  k = mean(x.^4) / max(mean(x.^2)^2, eps);
end

function [x, nrep] = despike(x)
% Replace impulsive samples (floe impacts): |x - movmedian| > 6 robust sigma.
    med = movmedian(x, 51);  r = x - med;
    s = 1.4826 * median(abs(r - median(r)));
    bad = abs(r) > 6*s;  x(bad) = med(bad);  nrep = nnz(bad);
end

function [psi, H] = estimate_heading(T, t_s, opts)
% psi(t) = psi0 + integral(yaw rate - offset) dt, anchored to the GNSS yaw proxy.
% Offset: burst mean (no net yaw inside a burst) + linear ramp if significant in
% 60 s block means. psi0 from fixes inside the burst, else the nearest fix.
    n = numel(t_s);  psi = nan(n, 1);
    H = struct('mode', "none", 'n_fix', 0, 'bias_dps', NaN);
    G = opts.gnss_yaw;
    if isempty(G) || ~ismember('Yaw', T.Properties.VariableNames), return; end
    if ismember('YawConfidence', G.Properties.VariableNames)
        G = G(G.YawConfidence >= opts.gnss_min_conf, :);
    end
    if isempty(G), return; end
    rate = opts.yaw_sign * opts.gyro_scale_dps * double(T.Yaw(:));  rate(isnan(rate)) = 0;
    dt = t_s(2) - t_s(1);
    b = mean(rate);
    L = round(60/dt);  nb = floor(n/L);
    if strcmpi(opts.yaw_bias_mode, 'per_burst_ramp') && nb >= 6
        rb = mean(reshape(rate(1:nb*L), L, nb), 1)';
        tb = ((0:nb-1)' + 0.5)*L*dt;  tc = tb - mean(tb);
        c = [ones(nb,1) tc] \ rb;
        s2 = sum((rb - [ones(nb,1) tc]*c).^2) / (nb - 2);
        if abs(c(2) / sqrt(s2/sum(tc.^2))) >= opts.ramp_t_min
            b = c(1) + c(2)*(t_s - mean(tb));
        end
    end
    gr = cumtrapz(rate - b) * dt;
    tf = seconds(notz(G.TimeUTC) - notz(T.DateTimeUTC(1)));
    yf = mod(opts.gnss_yaw_sense * G.YawEstimate_deg + opts.heading_offset_deg, 360);
    inb = tf >= 0 & tf <= t_s(end);
    H.n_fix = nnz(inb);
    if H.n_fix >= 1
        psi0 = circmean_deg(yf(inb) - interp1(t_s, gr, tf(inb)));
        H.mode = "gnss_anchor";
    else
        dist = max(-tf, tf - t_s(end));
        [gap, k] = min(dist);
        if isempty(gap) || gap > opts.gnss_max_gap_s, return; end
        if tf(k) < 0, psi0 = yf(k) - gr(1); else, psi0 = yf(k) - gr(end); end
        H.mode = "gnss_nearest";
    end
    psi = mod(psi0 + gr, 360);
    H.bias_dps = mean(b);
end

function t = notz(t)
    if isdatetime(t) && ~isempty(t.TimeZone), t.TimeZone = 'UTC'; t.TimeZone = ''; end
end

function m = circmean_deg(a)
    a = a(~isnan(a));
    if isempty(a), m = NaN; return; end
    m = mod(atan2d(mean(sind(a)), mean(cosd(a))), 360);
end

function s = fill_defaults(s, d)
    fn = fieldnames(d);
    for k = 1:numel(fn)
        if ~isfield(s, fn{k}), s.(fn{k}) = d.(fn{k}); end
    end
end

function v = iff(c, a, bfun)
    if c, v = a; else, v = bfun(); end
end