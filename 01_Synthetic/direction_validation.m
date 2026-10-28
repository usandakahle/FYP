%% direction_validation.m 
% direction check figures for the synthetic test 
% Builds the same synthetic packet as ocean_wave_spectrum.m + sharc_meas.m (same ERA5 hour,
% same RNG seed, same call order:  runs it through run_sharc_jack
% once, and plots ONLY the figures below. 
%    synthetic_dirfreq_<regime>   (a) direction per frequency, (b) spread per frequency, vs target
%    synthetic_qc_<regime>        (a) heave-slope coherence, (b) check ratio (both ~1 here)
%    synthetic_quadrature         surface elevation and slope, first 100 s
%    synthetic_target             (a) target E(f) and the realisation, (b) surface elevation
%                                   (same figure as ocean_wave_spectrum.m O1; MAKE_TARGET)
%    synthetic_dirspec_<regime>   target vs recovered E(f,theta) (only if MAKE_DIRSPEC = true)

clear; clc;
%% SETTINGS
REGIME       = "ice";            % "open" | "ice"
RNG_SEED     = 1;                % same seed as ocean_wave_spectrum.m -> same realisation
MAKE_DIRSPEC = false;            % V1 polar target vs SHARC (off: not needed, saves time)
MAKE_TARGET  = true;             %  target spectrum and realisation (synthetic_target.png)
ERA5_WAVE = 'Era5Wave.nc';
target_lat = -58.8;  target_lon = -0.69;
D        = 1800;
CI_FIXED = 1;
WCFG = struct('D', D, 'f_max', 0.6, 'ice_dist', 5000, 'a_coeff', 2.12e-3, 'b_coeff', 4.59e-2);
PCFG = struct('fs', 100, 'T_rec', D, 'scale_heave', 9.80665/8192, 'gyro_scale', 1/65.5, ...
    'n_acc_ug', 593, 'n_gyro_dps', 0.0352, 'gyro_bias', [0.3 -0.2 0.1]);
opts = struct('fs_imu','auto', 'gyro_scale_dps',PCFG.gyro_scale, 'n_acc_ug',PCFG.n_acc_ug, ...
    'n_gyro_dps',PCFG.n_gyro_dps, 'snr_min',3, 'noise_correct',true, 'f_min',0.03, 'f_max_band',0.35, ...
    'Nfft',512, 'M',32, 'tau',37, 'theta_offset_deg',0, 'yaw_sign',1, 'gnss_yaw_sense',1, ...
    'heading_offset_deg',0, 'jackknife',false, 'track',false, 'verbose',true, 'gnss_yaw',[]);
POLAR_FMAX = 0.30;

%% ERA5 sea state (peak-Hs hour, as ocean_wave_spectrum.m) 
rng(RNG_SEED);
FT = sharc_figtools();
REGIME = string(REGIME);
lat  = ncread(ERA5_WAVE, 'latitude');  lon = ncread(ERA5_WAVE, 'longitude');
time = ncread(ERA5_WAVE, 'valid_time');
swh  = ncread(ERA5_WAVE, 'swh');  pp1d = ncread(ERA5_WAVE, 'pp1d');
mwd  = ncread(ERA5_WAVE, 'mwd');  wdw  = ncread(ERA5_WAVE, 'wdw');
[~, ilat] = min(abs(lat - target_lat));  [~, ilon] = min(abs(lon - target_lon));
[~, itime] = max(squeeze(swh(ilon, ilat, :)));
target_time = time(itime);
Hs = swh(ilon, ilat, itime);  Tp = pp1d(ilon, ilat, itime);
mwd_deg = mwd(ilon, ilat, itime);  wdw_v = wdw(ilon, ilat, itime);

%% RUN
W = synth_truth(Hs, Tp, mwd_deg, wdw_v, CI_FIXED, WCFG);
t0 = datetime(1970,1,1, 'TimeZone','UTC') + seconds(double(target_time));
[Tpk, SIG] = synth_packet(W, REGIME, t0, PCFG);
[theta_sharc, Hm0_sharc, ~, ~, pd] = run_sharc_jack(Tpk, "synthetic " + REGIME, false, ...
    PCFG.scale_heave, 0, +1, [], opts);
R = synth_compare(W, REGIME, pd, opts.f_min);

dir_target = R.dir_target;  Czc = R.Czc;  chk_f = R.chk_f;
fr = R.f;  E_tgt = R.E_target;  dfp = pd.freq(2) - pd.freq(1);
t_imu = SIG.t;  eta_t = SIG.eta;  sE = SIG.sE;  sN = SIG.sN;
sigma = W.sigma;  s = W.s;  theta_mean = W.theta_mean;

rlab = "Open Water";  if REGIME == "ice", rlab = "In Ice"; end
reg = char(REGIME);
inb = pd.band(:);  fb = pd.freq(inb);
fsw = pd.freq(pd.sys_sw);  swr = [min(fsw) - dfp/2, max(fsw) + dfp/2];   % swell region (shaded)
cs  = [0.85 0.85 0.85];
fsl = FT.FIG.fs - 2;                                          % reference-line label size

%%  TARGET SPECTRUM AND REALISATION (as ocean_wave_spectrum.m)
if MAKE_TARGET
    fs2 = 2;  t2 = 0:1/fs2:D - 1/fs2;                         % realisation at 2 Hz, open water
    eta = zeros(size(t2));
    for i = 1:numel(W.f), eta = eta + W.a(i) * cos(2*pi*W.f(i)*t2 + W.alpha(i)); end
    [Pxx, f_est] = pwelch(eta, hamming(1024), 512, [], fs2);
    fig = FT.mkfig([1100 600]);
    tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    a1 = nexttile(tl);
    semilogy(a1, f_est, Pxx, '-', 'Color', FT.cb, 'LineWidth', FT.FIG.lw, 'DisplayName', 'Realisation (Welch)');
    hold(a1, 'on');  grid(a1, 'on');  box(a1, 'on');
    semilogy(a1, W.f, W.E_f, '-', 'Color', [0 0 0], 'LineWidth', 0.75, 'DisplayName', 'Target E(f)');
    yl = ylim(a1);  ylim(a1, [max(yl(1), 1e-6*max(W.E_f)) yl(2)]);
    xlim(a1, [0 0.5]);  xlabel(a1, 'Frequency (Hz)');  ylabel(a1, 'Spectral Density (m^2/Hz)');
    legend(a1, 'Location', 'northeast');
    title(a1, sprintf('Target (H_s %.2f m, T_p %.1f s, From %.0f^\\circ): (a) Spectrum, (b) Surface Elevation', ...
        Hs, Tp, mwd_deg), 'FontWeight', 'normal');
    a2 = nexttile(tl);  hold(a2, 'on');  grid(a2, 'on');  box(a2, 'on');
    iw_ = t2 <= 300;
    plot(a2, t2(iw_), eta(iw_), '-', 'Color', FT.cb, 'LineWidth', FT.FIG.lw);
    xlabel(a2, 'Time (s)');  ylabel(a2, 'Elevation \eta (m)');
    FT.save(fig, 'synthetic_target');
end

%%  TARGET vs RECOVERED
if MAKE_DIRSPEC
    tgc  = pd.theta_grid(:)';
    D_c  = max(cos(wrapToPi(deg2rad(90 - tgc) - theta_mean)/2), 0).^(2*s);
    E_tf = interp1(fr, E_tgt, pd.freq(:), 'linear', 0) .* D_c;
    fig = FT.mkfig([1100 480]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
    ax = nexttile(tl);
    FT.polar_E(ax, pd.freq(:), tgc(:), E_tf, POLAR_FMAX, ...
        struct('deg', dir_target, 'col', [0 0 0], 'ls', '--', 'name', 'Target'), FT.cmap);
    title(ax, '(a) Target', 'FontWeight', 'normal');
    ax = nexttile(tl);
    FT.polar_E(ax, pd.freq(:), tgc(:), Czc .* pd.D, POLAR_FMAX, ...
        struct('deg', {theta_sharc, dir_target}, 'col', {FT.cb, [0 0 0]}, 'ls', {'-', '--'}, ...
        'name', {'SHARC', 'Target'}), FT.cmap);
    title(ax, '(b) SHARC Pipeline', 'FontWeight', 'normal');
    FT.save(fig, ['synthetic_dirspec_' reg]);
end

%% DIRECTION AND SPREAD
th1 = dir_target + FT.wrap180(pd.theta1(:) - dir_target);     % unwrapped around the target
spr = rad2deg(sqrt(max(2*(1 - pd.r1(:)), 0)));
fig = FT.mkfig([1100 600]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
a1 = nexttile(tl);  hold(a1, 'on');  grid(a1, 'on');  box(a1, 'on');
shade(a1, swr, dir_target + [-90 90], cs, 'Swell Region');
plot(a1, fb, th1(inb), '-', 'Color', FT.cb, 'LineWidth', FT.FIG.lw, 'DisplayName', '\theta_1(f) (SHARC)');
yline(a1, dir_target, '--', 'Color', [0 0 0], 'DisplayName', sprintf('Target %.0f^\\circ', dir_target));
ylim(a1, dir_target + [-90 90]);  ylabel(a1, 'Mean Wave Direction (deg)');
legend(a1, 'Location', 'northeast');
title(a1, sprintf('Synthetic (%s): (a) Direction, (b) Spread per Frequency', rlab), 'FontWeight', 'normal');
a2 = nexttile(tl);  hold(a2, 'on');  grid(a2, 'on');  box(a2, 'on');
shade(a2, swr, [0 90], cs, '');
plot(a2, fb, spr(inb), '-', 'Color', FT.cb, 'LineWidth', FT.FIG.lw, 'DisplayName', '\sigma_1(f) (SHARC)');
yline(a2, rad2deg(sigma), '--', 'Color', [0 0 0], 'DisplayName', sprintf('Target %.0f^\\circ', rad2deg(sigma)));
ylim(a2, [0 90]);  ylabel(a2, 'Directional Spread (deg)');  xlabel(a2, 'Frequency (Hz)');
legend(a2, 'Location', 'northeast');
linkaxes([a1 a2], 'x');  xlim(a1, [0 opts.f_max_band]);  a1.XTickLabel = [];
FT.save(fig, ['synthetic_dirfreq_' reg]);

%% COHERENCE AND CHEK RATIO
fig = FT.mkfig([1100 600]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
a1 = nexttile(tl);  hold(a1, 'on');  grid(a1, 'on');  box(a1, 'on');
shade(a1, swr, [0 1.1], cs, '');
plot(a1, fb, pd.coh_f(inb), '-', 'Color', FT.cb, 'LineWidth', FT.FIG.lw);
yline(a1, 1, '--', 'Ideal (1)', 'Color', [0 0 0], 'LabelHorizontalAlignment', 'left', 'FontSize', fsl);
ylim(a1, [0 1.1]);  ylabel(a1, 'Coherence \gamma^2 (-)');
title(a1, sprintf(['Synthetic (%s): (a) Heave-Slope Coherence, ' ...
    '(b) Check Ratio (C_{xx}+C_{yy})/(k^2 C_{zz})'], rlab), 'FontWeight', 'normal');
a2 = nexttile(tl);  hold(a2, 'on');  grid(a2, 'on');  box(a2, 'on');
shade(a2, swr, [0 2], cs, '');
plot(a2, fb, chk_f(inb), '-', 'Color', FT.cb, 'LineWidth', FT.FIG.lw);
yline(a2, 1, '--', 'Ideal (1)', 'Color', [0 0 0], 'LabelHorizontalAlignment', 'left', 'FontSize', fsl);
ylim(a2, [0 2]);  ylabel(a2, 'Check Ratio (-)');  xlabel(a2, 'Frequency (Hz)');
linkaxes([a1 a2], 'x');  xlim(a1, [0 opts.f_max_band]);  a1.XTickLabel = [];
FT.save(fig, ['synthetic_qc_' reg]);

%% SURFACE ELEVATION AND SLOPE
iw = t_imu <= 100;
fig = FT.mkfig([1100 550]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
a1 = nexttile(tl);  hold(a1, 'on');  grid(a1, 'on');  box(a1, 'on');
plot(a1, t_imu(iw), eta_t(iw), '-', 'Color', FT.cb, 'LineWidth', FT.FIG.lw);
ylabel(a1, 'Elevation \eta (m)');  a1.XTickLabel = [];
title(a1, sprintf('Synthetic (%s): (a) Surface Elevation, (b) Surface Slope', rlab), 'FontWeight', 'normal');
a2 = nexttile(tl);  hold(a2, 'on');  grid(a2, 'on');  box(a2, 'on');
plot(a2, t_imu(iw), 1e3*sN(iw), '-', 'Color', FT.cb, 'LineWidth', FT.FIG.lw, 'DisplayName', 's_N (North)');
plot(a2, t_imu(iw), 1e3*sE(iw), '-', 'Color', FT.c2, 'LineWidth', FT.FIG.lw, 'DisplayName', 's_E (East)');
ylabel(a2, 'Slope (mrad)');  xlabel(a2, 'Time (s)');  legend(a2, 'Location', 'northeast');
linkaxes([a1 a2], 'x');  xlim(a1, [0 100]);
FT.save(fig, 'synthetic_quadrature');

fprintf('Figures in %s/ (synthetic_target, synthetic_dirfreq_%s, synthetic_qc_%s, synthetic_quadrature)\n', ...
    FT.FIG.outdir, reg, reg);

%% LOCAL FUNCTION
function shade(ax, xr, yr, c, name)
% Light band marking the swell region the direction is taken from (in the legend if named).
    h = patch(ax, [xr fliplr(xr)], [yr(1) yr(1) yr(2) yr(2)], c, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
    if strlength(name) > 0, h.DisplayName = name; else, h.HandleVisibility = 'off'; end
end
