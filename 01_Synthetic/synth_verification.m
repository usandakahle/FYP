%% synth_verification.m 
% synthetic verification table (E1, E2a-c, E3)
% Runs N_SEEDS realisations per regime through the FINAL pipeline
% (run_sharc_jack)
% and writes one row per seed to 03TestFigures/synthetic/synth_verification.csv,

% Figures:
%   verif_packet_<case>.png       packet figure (field form) for seed FIG_SEED of each case
%   synth_verification_errors.png (a) Hm0 error, (b) direction error vs realised, per seed

% Tests
%   E1  Hm0 vs the REALISED and the TARGET value over the same band      
%   E2a clean cross-spectra: truth eta, sE, sN (no noise, no filters), same segments and bins
%       as the pipeline. Direction vs the realised direction 
%   E2b fusion: fused slope (pd.ts, needs run_sharc_jack option return_ts) vs true slope,
%       RMS error in the 0.03-0.5 Hz band (requirement: < 0.057 deg)
%   E2c full chain: direction and spread vs REALISED (verification) and vs TARGET 
% , with sigma_theta and the jackknife se
%
% Reference choice: each component gets ONE random direction (Monte Carlo sample of D(theta)),
% so a realisation has its own mean direction. The chain is verified against that realised
% direction; target - realised is sampling scatter, not chain error.

clear; clc;
FT = sharc_figtools();

%% SETTINGS
N_SEEDS   = 20;
REGIMES   = ["open" "ice"];
RUN_E3    = true;  E3_REGIME = "open";
PIPE      = @run_sharc_jack;          % final pipeline
ERA5_WAVE = 'Era5Wave.nc';
target_lat = -58.8;  target_lon = -0.69;

% Fallback sea state if the ERA5 file is not on the path (values of the report figures)
FALLBACK  = struct('Hs', 5.12, 'Tp', 12.7, 'mwd', 311, 'wdw', deg2rad(26));
CI_ICE    = 1;                        % ice concentration for the "ice" regime
WCFG = struct('D', 1800, 'f_max', 0.6, 'ice_dist', 5000, 'a_coeff', 2.12e-3, 'b_coeff', 4.59e-2);
PCFG = struct('fs', 100, 'T_rec', WCFG.D, 'scale_heave', 9.80665/8192, 'gyro_scale', 1/65.5,'n_acc_ug',...
    593, 'n_gyro_dps', 0.0352, 'gyro_bias', [0.3 -0.2 0.1]);

% E3 wind sea, added to the ERA5 swell
WS = struct('Hs', 2.0, 'Tp', 6.0, 'dtheta', 90, 'wdw', deg2rad(30));

% final pipeline options (same as batchBias)
opts = struct('fs_imu','auto', 'gyro_scale_dps',PCFG.gyro_scale, 'n_acc_ug',PCFG.n_acc_ug, ...
    'n_gyro_dps',PCFG.n_gyro_dps, 'snr_min',3, 'noise_correct',true, 'f_min',0.03, 'f_max_band',0.35, ...
    'Nfft',512, 'M',32, 'tau',37, 'theta_offset_deg',0, 'yaw_sign',1, 'gnss_yaw_sense',1, ...
    'heading_offset_deg',0, 'jackknife',true, 'track',false, 'verbose',false, 'gnss_yaw',[], ...
    'seg_mode',"robust", 'kf_class_order',"impacts_first", 'hm0_range',"band", ...
    'dir_fmax_wind',false, 'return_ts',true);

% dir_fmax_wind = false: the synthetic sea has no ERA5 wind; the fixed 0.15 Hz cut off is used.
F_BAND = [0.03 0.5];          % fusion check band (Hz)
CI_P   = 0.90;                % E(f) confidence level
FIG_SEED   = 1;               % seed whose packet figure is saved for each case (0: none)
POLAR_FMAX = 0.30;            % Hz, outer ring of the polar plot

%% Sea state 
if isfile(ERA5_WAVE)
    lat = ncread(ERA5_WAVE,'latitude');  lon = ncread(ERA5_WAVE,'longitude');
    swh = ncread(ERA5_WAVE,'swh');  pp1d = ncread(ERA5_WAVE,'pp1d');
    mwd = ncread(ERA5_WAVE,'mwd');  wdw = ncread(ERA5_WAVE,'wdw');
    [~, ilat] = min(abs(lat - target_lat));  [~, ilon] = min(abs(lon - target_lon));
    [~, it] = max(squeeze(swh(ilon, ilat, :)));
    SS = struct('Hs', double(swh(ilon,ilat,it)), 'Tp', double(pp1d(ilon,ilat,it)), ...
                'mwd', double(mwd(ilon,ilat,it)), 'wdw', double(wdw(ilon,ilat,it)));
else
    warning('%s not found: using the fallback sea state.', ERA5_WAVE);  SS = FALLBACK;
end
fprintf('Sea state: Hs %.2f m, Tp %.1f s, from %.0f deg, wdw %.3f rad (%.0f deg)\n', ...
    SS.Hs, SS.Tp, SS.mwd, SS.wdw, rad2deg(SS.wdw));

%% Run
rows = {};
cases = [repmat(REGIMES(:), 1, 1)];
if RUN_E3, cases = [cases; "two_" + E3_REGIME]; end
t0 = datetime(2022,7,21,12,0,0, 'TimeZone','UTC');
tic;
for c = 1:numel(cases)
    cs = cases(c);  two = startsWith(cs, "two_");  reg = erase(cs, "two_");
    ci = 0;  if reg == "ice", ci = CI_ICE; end
    for sd = 1:N_SEEDS
        rng(sd);
        W1 = synth_truth(SS.Hs, SS.Tp, SS.mwd, SS.wdw, ci, WCFG);
        W  = W1;  nsw = numel(W1.f);
        if two
            W2 = synth_truth(WS.Hs, WS.Tp, mod(SS.mwd + WS.dtheta, 360), WS.wdw, ci, WCFG);
            W  = merge_fields(W1, W2);
        end
        [Tpk, SIG] = synth_packet(W, reg, t0, PCFG);
        try
            [th, Hm0, Tm01, ~, pd] = PIPE(Tpk, "verif " + cs + " seed " + sd, false, ...
                PCFG.scale_heave, 0, +1, [], opts);
        catch ME
            fprintf(2, '%s seed %d failed: %s\n', cs, sd, ME.message);  continue
        end
        r = struct('case', cs, 'seed', sd);
        r = e1(r, W, reg, pd, Hm0, CI_P, two, W1, nsw);
        r = e2a(r, SIG, W, reg, pd, opts, nsw);
        r = e2b(r, SIG, pd, F_BAND);
        r = e2c(r, W, reg, pd, th, nsw);
        rows{end+1} = r; %#ok<SAGROW>
        if sd == FIG_SEED
            try
                verif_packet_fig(FT, pd, cs, sd, Hm0, th, r, W, reg, POLAR_FMAX);
            catch ME
                warning('packet figure %s seed %d not made: %s', cs, sd, ME.message);
            end
        end
        fprintf('[%-8s %2d] Hm0 %.2f (real %.2f)  th %.1f (real %.1f, tgt %.1f)  clean %+.2f  fuse %.3f deg  %.0f s\n', ...
            cs, sd, Hm0, r.Hm0_real, th, r.dir_real, r.dir_target, r.e2a_err_deg, r.e2b_rms_band_deg, toc);
    end
end
R = struct2table([rows{:}]);

%% CSV 
outd = FT.FIG.outdir;  if ~exist(outd, 'dir'), mkdir(outd); end
writetable(R, fullfile(outd, 'synth_verification.csv'));

%% Summary
cm  = @(e) rad2deg(angle(mean(exp(1i*deg2rad(e(isfinite(e)))))));
cs_ = @(e) rad2deg(sqrt(-2*log(abs(mean(exp(1i*deg2rad(e(isfinite(e)))))))));
S = table();
for c = 1:numel(cases)
    k = R.case == cases(c);  X = R(k, :);  n = height(X);
    s = struct();
    s.case = cases(c);  s.n = n;  s.n_valid = nnz(X.valid == 1);
    s.Hm0_err_real_mean_pct = mean(X.Hm0_err_real_pct);
    s.Hm0_err_real_se_pct   = std(X.Hm0_err_real_pct)/sqrt(n);
    s.Hm0_err_real_maxabs_pct = max(abs(X.Hm0_err_real_pct));
    s.Hm0_err_tgt_mean_pct  = mean(X.Hm0_err_tgt_pct);
    s.Hm0_err_tgt_std_pct   = std(X.Hm0_err_tgt_pct);
    s.Hm0_err_tgt_maxabs_pct = max(abs(X.Hm0_err_tgt_pct));
    s.Ef_in_ci_mean_pct     = mean(X.Ef_in_ci_pct);
    s.Ef_in_ci_min_pct      = min(X.Ef_in_ci_pct);
    s.n_bimodal             = nnz(X.bimodal == 1);
    s.e2a_maxabs_deg        = max(abs(X.e2a_err_deg));
    s.e2b_rms_band_mean_deg = mean(X.e2b_rms_band_deg);
    s.e2b_rms_band_max_deg  = max(X.e2b_rms_band_deg);
    s.e2b_slope_rms_deg     = mean(X.e2b_slope_rms_deg);
    s.dir_err_real_mean_deg = cm(X.dir_err_real);
    s.dir_err_real_cstd_deg = cs_(X.dir_err_real);
    s.dir_err_real_maxabs   = max(abs(X.dir_err_real));
    s.dir_err_tgt_mean_deg  = cm(X.dir_err_tgt);
    s.dir_err_tgt_cstd_deg  = cs_(X.dir_err_tgt);
    s.sampling_cstd_deg     = cs_(FT.wrap180(X.dir_real - X.dir_target));
    s.sigma_theta_mean_deg  = mean(X.sigma_theta, 'omitnan');
    s.jk_se_mean_deg        = mean(X.jk_se, 'omitnan');
    s.crit_dir_deg          = 2*s.sigma_theta_mean_deg/sqrt(n);
    s.spread_err_real_mean  = mean(X.spread_err_real, 'omitnan');
    s.spread_err_tgt_mean   = mean(X.spread_err_tgt, 'omitnan');
    s.r1_mean               = mean(X.r1_sys, 'omitnan');
    s.n_bins_sys_mean       = mean(X.n_bins_sys, 'omitnan');
    s.band_dir_err_mean_deg = cm(X.band_err_real);
    S = [S; struct2table(s)]; %#ok<AGROW>
end
writetable(S, fullfile(outd, 'synth_verification_summary.csv'));
disp(S);
fprintf('Rows in %s, summary in synth_verification_summary.csv\n', fullfile(outd, 'synth_verification.csv'));

%%  FIGURE: per seed errors
cols = {FT.cb, FT.c2, [0.20 0.60 0.30]};
fig = FT.mkfig([1100 600]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
a1 = nexttile(tl);  hold(a1, 'on');  grid(a1, 'on');  box(a1, 'on');
a2 = nexttile(tl);  hold(a2, 'on');  grid(a2, 'on');  box(a2, 'on');
for c = 1:numel(cases)
    k = R.case == cases(c);  col = cols{1 + mod(c - 1, numel(cols))};
    plot(a1, R.seed(k), R.Hm0_err_real_pct(k), '-o', 'Color', col, 'LineWidth', FT.FIG.lw, ...
        'MarkerSize', 4, 'MarkerFaceColor', col, 'DisplayName', case_label(cases(c)));
    plot(a2, R.seed(k), R.dir_err_real(k), '-o', 'Color', col, 'LineWidth', FT.FIG.lw, ...
        'MarkerSize', 4, 'MarkerFaceColor', col, 'DisplayName', case_label(cases(c)));
end
yline(a1, 0, '--', 'Color', FT.cg, 'HandleVisibility', 'off');
yline(a2, 0, '--', 'Color', FT.cg, 'HandleVisibility', 'off');
ylabel(a1, 'H_{m0} Error vs Realised (%)');  a1.XTickLabel = [];
title(a1, 'Synthetic Verification: (a) H_{m0} Error, (b) Direction Error', 'FontWeight', 'normal');
legend(a1, 'Location', 'best');
ylabel(a2, 'Direction Error vs Realised (deg)');  xlabel(a2, 'Seed (-)');
linkaxes([a1 a2], 'x');  xlim(a1, [0.5 N_SEEDS + 0.5]);
FT.save(fig, 'synth_verification_errors');
fprintf('Figures in %s/ (verif_packet_*.png, synth_verification_errors.png)\n', outd);

%% ======================= local functions =======================
function W = merge_fields(W1, W2)
% Concatenate the components of two wave systems (swell first). synth_packet sums over
% components, so duplicated frequencies are fine. Spectra are kept on W1's grid.
    W = W1;
    for fn = ["f" "a" "a_ice" "alpha" "theta_i"]
        W.(fn) = [W1.(fn) W2.(fn)];
    end
    W.E_f = W1.E_f + W2.E_f;  W.Ef_ice = W1.Ef_ice + W2.Ef_ice;  W.f_grid = W1.f;
end

function [amp, E_t, fgrid] = pick(W, reg)
    amp = W.a(:);  E_t = W.E_f(:);
    if reg == "ice", amp = W.a_ice(:);  E_t = W.Ef_ice(:); end
    fgrid = W.f(:);  if isfield(W, 'f_grid'), fgrid = W.f_grid(:); end
end

function r = e1(r, W, reg, pd, Hm0, p, two, W1, nsw) %#ok<INUSD>
% Hm0 over the pipeline's own band; E(f) inside the chi-square interval; false second peak.
    [amp, E_t, fg] = pick(W, reg);  fr = W.f(:);
    dfp = pd.freq(2) - pd.freq(1);
    fb  = pd.freq(pd.band);
    r.f_lo = NaN; r.f_hi = NaN; r.Hm0 = Hm0; r.Hm0_real = NaN; r.Hm0_target = NaN;
    r.Hm0_err_real_pct = NaN; r.Hm0_err_tgt_pct = NaN; r.Ef_in_ci_pct = NaN;
    if ~isempty(fb)
        lo = fb(1) - dfp/2;  hi = fb(end) + dfp/2;  r.f_lo = fb(1);  r.f_hi = fb(end);
        m = fr >= lo & fr <= hi;   mg = fg >= lo & fg <= hi;
        r.Hm0_real   = 4*sqrt(sum(amp(m).^2/2));
        r.Hm0_target = 4*sqrt(trapz(fg(mg), E_t(mg)));
        r.Hm0_err_real_pct = 100*(Hm0/r.Hm0_real - 1);
        r.Hm0_err_tgt_pct  = 100*(Hm0/r.Hm0_target - 1);
        % E(f): bin-averaged target vs noise-corrected estimate, nu = 2 x segments used
        Czc = max(pd.Czz(:) - pd.Nzz(:), 0);
        nu = 2*max(pd.nseg_used, 1);
        if isfield(pd, 'robust_neff') && ~isempty(pd.robust_neff) && any(isfinite(pd.robust_neff(:)))
            nu = 2*mean(pd.robust_neff(:), 'omitnan');               % robust mode: effective segments
        end
        lim = [nu/chi2inv((1+p)/2, nu), nu/chi2inv((1-p)/2, nu)];
        ib = find(pd.band);  inside = false(size(ib));
        for j = 1:numel(ib)
            fj = pd.freq(ib(j));  mj = fg >= fj - dfp/2 & fg < fj + dfp/2;
            Et = mean(E_t(mj));
            inside(j) = Czc(ib(j)) >= lim(1)*Et && Czc(ib(j)) <= lim(2)*Et;
        end
        r.Ef_in_ci_pct = 100*mean(inside);
    end
    r.bimodal = double(isfield(pd, 'bimodal') && any(pd.bimodal));
    r.valid = double(pd.dir_valid);
    r.invalid_reason = string(pd.invalid_reason);
    r.nseg_used = pd.nseg_used;
end

function [dir, r1] = realised_dir(amp, th_i, fr, msk)
% Energy-weighted realised direction FROM (compass) and r1 of the components in msk.
    z = sum(amp(msk).^2 .* exp(1i*th_i(msk)));
    dir = mod(90 - rad2deg(angle(z)), 360);
    r1 = abs(z) / max(sum(amp(msk).^2), eps);
end

function msk = comp_mask(fr, pd, dfp, which)
    msk = false(size(fr));
    if ~any(pd.(which)), return; end
    f_ = pd.freq(pd.(which));
    msk = fr >= min(f_) - dfp/2 & fr <= max(f_) + dfp/2;
end

function r = e2a(r, SIG, W, reg, pd, opts, nsw)
% Clean cross-spectra: truth signals, plain downsampling (content < 0.6 Hz, no aliasing at
% 3.125 Hz), same Hann/Bartlett segments and bins as the pipeline, no noise, no filters.
    M = opts.M;  N = opts.Nfft;
    z  = SIG.eta(1:M:end);  sE = SIG.sE(1:M:end);  sN = SIG.sN(1:M:end);
    K = floor(numel(z)/N);  w = hann(N, 'periodic')';
    QE = zeros(1, N/2+1);  QN = QE;
    for k = 1:K
        ii = (k-1)*N + (1:N);
        Z = fft(z(ii).*w);  E = fft(sE(ii).*w);  Nn = fft(sN(ii).*w);
        QE = QE + imag(conj(Z(1:N/2+1)) .* E(1:N/2+1));
        QN = QN + imag(conj(Z(1:N/2+1)) .* Nn(1:N/2+1));
    end
    r.e2a_err_deg = NaN;  r.e2a_dir = NaN;
    if any(pd.sys_sw)
        fc = (0:N/2)*pd.fs/M/N;                                     % clean-spectrum bins
        m = ismembertol(fc, pd.freq(pd.sys_sw), 1e-6);              % same bins as the pipeline
        r.e2a_dir = mod(90 - atan2d(sum(QN(m)), sum(QE(m))), 360);
        [amp, ~] = pick(W, reg);  fr = W.f(:);  dfp = pd.freq(2) - pd.freq(1);
        msk = comp_mask(fr, pd, dfp, 'sys_sw');
        dr = realised_dir(amp, W.theta_i(:), fr, msk);
        r.e2a_err_deg = mod(r.e2a_dir - dr + 180, 360) - 180;
    end
end

function r = e2b(r, SIG, pd, fb)
% Fusion: fused slope (gyro convention Sx = pitch ~ s_N, Sy = roll ~ s_E) vs truth, in band.
    r.e2b_rms_band_deg = NaN;  r.e2b_slope_rms_deg = NaN;  r.e2b_rms_total_deg = NaN;
    if ~isfield(pd, 'ts'), return; end
    n = min(numel(pd.ts.Sx), numel(SIG.sN));
    eN = pd.ts.Sx(1:n) - SIG.sN(1:n)';  eE = pd.ts.Sy(1:n) - SIG.sE(1:n)';
    fs = 1/median(diff(pd.ts.t));
    bp = @(x) band(x - mean(x), fs, fb);
    r.e2b_rms_band_deg  = rad2deg(sqrt(mean(bp(eN).^2 + bp(eE).^2)/2));
    r.e2b_rms_total_deg = rad2deg(sqrt(mean((eN - mean(eN)).^2 + (eE - mean(eE)).^2)/2));
    r.e2b_slope_rms_deg = rad2deg(sqrt(mean(bp(SIG.sN(1:n)').^2 + bp(SIG.sE(1:n)').^2)/2));
end

function y = band(x, fs, fb)
% Ideal band-pass by FFT mask (whole record).
    x = x(:);  n = numel(x);  X = fft(x);
    f = (0:n-1)'*fs/n;  f = min(f, fs - f);
    X(f < fb(1) | f > fb(2)) = 0;
    y = real(ifft(X));
end

function r = e2c(r, W, reg, pd, th, nsw)
% Full chain direction and spread vs realised (swell components only) and vs target.
    [amp, ~] = pick(W, reg);  fr = W.f(:);  th_i = W.theta_i(:);
    dfp = pd.freq(2) - pd.freq(1);
    sw = true(size(fr));  sw(nsw+1:end) = false;               % swell components (E3)
    msk = comp_mask(fr, pd, dfp, 'sys_sw') & sw;
    r.theta = th;  r.dir_target = mod(90 - rad2deg(W.theta_mean), 360);
    [r.dir_real, r1r] = realised_dir(amp, th_i, fr, msk);
    r.dir_err_real = mod(th - r.dir_real + 180, 360) - 180;
    r.dir_err_tgt  = mod(th - r.dir_target + 180, 360) - 180;
    r.spread_target = rad2deg(W.sigma);
    r.spread_real   = rad2deg(sqrt(max(2*(1 - r1r), 0)));
    r.spread = pd.spread_sys;
    r.spread_err_real = pd.spread_sys - r.spread_real;
    r.spread_err_tgt  = pd.spread_sys - r.spread_target;
    r.r1_sys = pd.r1_sys;  r.n_bins_sys = pd.n_bins_sys;
    r.sigma_theta = pd.sigma_theta;  r.jk_se = pd.jk_se;
    % band-mean direction (no partition) vs realised swell: shows what the partition buys (E3)
    r.band_err_real = mod(pd.theta_band - r.dir_real + 180, 360) - 180;
end

function verif_packet_fig(FT, pd, cs, sd, Hm0, th, r, W, reg, fmax)
% Packet figure in the field form (sharc_figtools packet_fig) for one verification run.
    [~, E_t, fg] = pick(W, reg);
    head = {sprintf('Synthetic Verification: %s, Seed %d', case_label(cs), sd), ...
        sprintf(['H_{m0} = %.2f m,   \\theta_{alg} = %.0f^\\circ,   Realised = %.0f^\\circ,   ' ...
        'Target = %.0f^\\circ,   r_1 = %.2f'], Hm0, th, r.dir_real, r.dir_target, pd.r1_sys)};
    dirs = struct('deg', {th, r.dir_real, r.dir_target}, 'col', {FT.cb, FT.c2, [0 0 0]}, ...
        'ls', {'-', '-', '--'}, 'name', {'SHARC', 'Realised', 'Target'});
    fig = FT.packet_fig(pd, head, dirs, fmax, struct('f', fg, 'E', E_t, 'name', 'Target E(f)'));
    FT.save(fig, char("verif_packet_" + cs));
    close(fig);
end

function s = case_label(cs)
    switch string(cs)
        case "open",     s = 'Open Water';
        case "ice",      s = 'In Ice';
        case "two_open", s = 'Two Systems (Open Water)';
        case "two_ice",  s = 'Two Systems (In Ice)';
        otherwise,       s = char(cs);
    end
end
