%% batchBias.m  
% swell direction and the storm diagnosis, as report figures

% Runs the pipeline on every packet of one buoy, tracks the swell direction across packets
% (Kalman filter + RTS smoother, each packet weighted by its own uncertainty) and, after a
% storm onset fixed by rule, adds the heading proxy psi as a second, gated measurement.
% Reference: LP8 only. 
%   F1 <buoy>_track_proxy      LP8, psi, tracked theta (SHARC body), tracked theta with heading (psi after onset)
%   F2 <buoy>_onset_rule       (a) max |yaw rate|, (b) kurtosis > 2 Hz, with thresholds
%   F5 <buoy>_coherence        (a) heave slope coherence over the swell region, (b) |theta - LP8|
%   CSV <buoy>_results.csv     one row per packet (Hm0, Tm01, directions, errors, wind, ice, coherence)
%   F3 <buoy>_packets/<buoy>_packet_<ddHHMM>  for EVERY packet
%   F4 <buoy>_wind/<buoy>_wind_<ddHHMM>  for the same packets: (a) U10 over the deployment with a
%                              green line at the packet time, (b) the packet's wind direction on a
%                              compass (N, E, S, W): line from the centre towards where the wind
%                              comes FROM, length = U10

% SHARC quantities in the buoy colour (SB04 blue, SB06 green) 
% psi grey; LP8 black dashed; wind black. 
% The tracked series is a FILTERED estimate: in the storm it interpolates between valid
% packets and, with psi, is led mainly by psi. 

clc;
%% SETTINGS
BUOY_ID        = "SB04";                                   
PACKET_DIR     = BUOY_ID + "/" + lower(BUOY_ID) + "_packets";
GNSS_POS_FILE  = BUOY_ID + "/" + BUOY_ID + "GNSSTable.csv";
YAW_STATE_FILE = BUOY_ID + "_yaw_state.mat";
REF_XLSX       = "reference_directions.xlsx";              % LP8 Mdir (sheet 'reference')
ERA5_NC        = "Era5WindIce.nc";                         % u10, v10 (wind for F4)
REF_DIR_WIN_MIN = 60;
windowStart    = datetime(2022,7,1, 'TimeZone','UTC');
windowEnd      = datetime(2022,8,1, 'TimeZone','UTC');
SCALE_HEAVE    = 9.80665 / 8192;
SIGN_QUAD      = +1;
PIPELINE       = @run_sharc_jack;   % pipeline used 
% tracker
Q_RATE      = 10;      % deg / sqrt(h): swell direction random walk (LP8 rate statistics)
SIG_FLOOR   = 10;      % deg, smallest measurement std for a valid packet (bias is not in sigma_theta)
SIG_INVALID = 60;      % deg, measurement std for a packet the pipeline marks invalid
GATE        = 3;       % innovations beyond GATE sigma are rejected
GATE_RESET  = 3;       % this many consistent rejections in a row restart the track
% heading proxy after the onset
ONSET_MODE  = "yaw";   % "yaw": max |yaw rate| 
YAW_ONSET   = 5;       % deg/s
KURT_ONSET  = 3.5;
KURT_GAUSS  = 3;       % reference line: Gaussian motion
ONSET_RUN   = 3;       % the condition must hold for this many packets in a row
ONSET_TIME  = datetime(2022,7,22,16,52,0);   % used when ONSET_MODE = "time" (UTC)
F2_DROP_START = 1;     % packets left out of the onset figure at the start (deployment transient)
F2_DROP_END   = 1;     % ... and at the end (recovery transient); set 0 to show them (plot only)
R1_VALID      = 0.85;  % pipeline r1 validity (reference line on the coherence figure, F5)
PSI_SIGMA   = 30;      % deg, std given to the heading proxy as a wave direction measurement
PSI_OFFSET  = 0;       % deg, added to psi (0: psi used as a FROM direction)
% packet figures (F3): "all" = one figure per packet (in <outdir>/<buoy>_packets);

PACKET_FIGS      = "all";
PACKET_FIG_TIMES = datetime.empty;          % e.g. [datetime(2022,7,20,15,59,0) datetime(2022,7,24,5,54,0)]
POLAR_FMAX  = 0.30;    % Hz, outer ring of the polar plot
GREEN       = [0.20 0.60 0.30];   % packet marker on the wind figures


opts = struct('fs_imu','auto', 'gyro_scale_dps',1/65.5, 'n_acc_ug',593, 'n_gyro_dps',0.0352, ...
    'snr_min',3, 'noise_correct',true, 'f_min',0.03, 'f_max_band',0.35, 'Nfft',512, 'M',32, ...
    'tau',37, 'theta_offset_deg',0, 'yaw_sign',1, 'gnss_yaw_sense',1, 'heading_offset_deg',0, ...
    'jackknife',false, 'track',false, 'verbose',false);

%% STYLE                       
FIG = struct('outdir', "OutputFigures", 'dpi', 200, 'lw', 1.0, 'ms', 3, 'fs', 16);
% fs: base font size (pt) for ticks, labels, titles and text. The figures are ~1100 px wide and
% are scaled to ~0.55 of that on an A4 page, so 16 pt prints at roughly 9 pt.
set_plot_defaults(FIG.fs);
tag = lower(char(BUOY_ID));
cb  = [0.00 0.45 0.74];  if upper(BUOY_ID) == "SB06", cb = [0.20 0.60 0.30]; end   % buoy colour
c2  = [0.85 0.33 0.10];  cg = [0.50 0.50 0.50];
cband = [0.30 0.80 0.30];  if upper(BUOY_ID) == "SB06", cband = [0.30 0.55 0.90]; end % SNR band shade
cmap = parula(256);

%% INPUTS 
G = [];
if isfile(YAW_STATE_FILE)
    S_ = load(YAW_STATE_FILE);
    G = table(S_.Results.TimeUTC, S_.Results.YawEstimate_deg, 'VariableNames', {'TimeUTC','YawEstimate_deg'});
    if isempty(G.TimeUTC.TimeZone), G.TimeUTC.TimeZone = 'UTC'; end
    G = sortrows(G(~isnan(G.YawEstimate_deg) & ~isnat(G.TimeUTC), :), 'TimeUTC');
else
    warning('%s not found: body frame only.', YAW_STATE_FILE);
end
opts.gnss_yaw = G;
GPOS = [];  if isfile(GNSS_POS_FILE), GPOS = readtable(GNSS_POS_FILE); end
R = [];
if isfile(REF_XLSX)
    R = readtable(REF_XLSX, 'Sheet', 'reference', 'TextType', 'string');
    if ~isdatetime(R.TimeUTC), R.TimeUTC = datetime(R.TimeUTC); end
    R.TimeUTC.TimeZone = '';
end
E5 = era5_load(ERA5_NC, {'u10n','v10n','siconc'});        % neutral 10 m wind, sea-ice concentration
UV = {'u10n','v10n'};
if isempty(E5) || ~isfield(E5, 'u10n')                     % fall back to u10/v10 if that is what the file has
    E5 = era5_load(ERA5_NC, {'u10','v10','siconc'});  UV = {'u10','v10'};
end

files = dir(fullfile(PACKET_DIR, "packet_*.csv"));
t0s = NaT(numel(files), 1, 'TimeZone', 'UTC');
for k = 1:numel(files)
    tok = regexp(files(k).name, 'packet_(\d+)_(.+)\.csv', 'tokens', 'once');
    if ~isempty(tok), t0s(k) = datetime(tok{2}, 'InputFormat','yyyy-MM-dd_HH-mm-ss', 'TimeZone','UTC'); end
end
keep = t0s >= windowStart & t0s <= windowEnd;
files = files(keep);  [~, o] = sort(t0s(keep));  files = files(o);
nP = numel(files);
fprintf('batchBias %s: %d packets\n', BUOY_ID, nP);

%% RUN
vn = {'theta','sig','valid','kurt','ref_lp8','psi','yaw_max','Hm0','r1','f_sw','U10','wind_from','Tm01','ice','coh'};
V = nan(nP, numel(vn));  tm = NaT(nP, 1);  SP = cell(nP, 1);
tic;
for f = 1:nP
    try
        Tw = readtable(fullfile(files(f).folder, files(f).name));
        if ~isdatetime(Tw.DateTimeUTC)
            Tw.DateTimeUTC = datetime(Tw.DateTimeUTC, 'InputFormat','dd-MMM-yyyy HH:mm:ss', 'TimeZone','UTC');
        elseif isempty(Tw.DateTimeUTC.TimeZone)
            Tw.DateTimeUTC.TimeZone = 'UTC';
        end
        [~, ~, ~, ~, p] = PIPELINE(Tw, files(f).name, false, SCALE_HEAVE, 0, SIGN_QUAD, [], opts);
        tm(f) = p.t_mid;
        [la, lo] = pos_at(GPOS, p.t_mid);
        u = era5_at(E5, UV{1}, p.t_mid, la, lo);  w = era5_at(E5, UV{2}, p.t_mid, la, lo);
        V(f,:) = [p.theta, p.sigma_theta, p.dir_valid, p.kurt_hf, lp8_ref(R, p.t_mid, REF_DIR_WIN_MIN), ...
                  p.psi_mid, p.yaw_max, p.Hm0, p.r1_sys, p.f_sw, hypot(u, w), mod(atan2d(-u, -w), 360), ...
                  p.Tm01, era5_at(E5, 'siconc', p.t_mid, la, lo), p.coh_sys];
        SP{f} = struct('freq', p.freq, 'Czz', p.Czz, 'Nzz', p.Nzz, 'band', p.band, 'D', p.D, ...
            'tg', p.theta_grid, 'theta', p.theta, 'Hm0', p.Hm0, 'r1', p.r1_sys, 'valid', p.dir_valid, ...
            'dir_fmax', p.dir_fmax, 't', p.t_mid);
    catch ME
        fprintf(2, '%s failed: %s\n', files(f).name, ME.message);
    end
    if mod(f, 10) == 0 || f == nP, fprintf('[%3d/%3d] %.0f s\n', f, nP, toc); end
end
C = array2table(V, 'VariableNames', vn);  C.t = tm;
ok_ = ~isnat(C.t);  C = C(ok_, :);  SP = SP(ok_);
if height(C) == 0
    error(['batchBias: no packets to plot. %d packet files found in %s between %s and %s; ' ...
        'if 0, run from the folder that contains %s/ or fix PACKET_DIR; ' ...
        'otherwise every packet failed in the pipeline (see the red "failed" lines above).'], ...
        nP, PACKET_DIR, string(windowStart), string(windowEnd), BUOY_ID);
end
v1 = C.valid == 1;
x  = tday(C.t);                                     % day of July 2022 (UTC)

%% TRACKER 
Rm = SIG_INVALID * ones(height(C), 1);
Rm(v1) = max(C.sig(v1), SIG_FLOOR);
Rm(~isfinite(Rm)) = SIG_INVALID;
[C.track, C.track_sd, C.gated] = track_dir(C.t, C.theta, Rm, Q_RATE, GATE, GATE_RESET);
switch ONSET_MODE
    case "yaw",  hit = C.yaw_max > YAW_ONSET;
    case "kurt", hit = C.kurt > KURT_ONSET;
    otherwise,   hit = false(height(C), 1);
end
k_on = NaN;
if ONSET_MODE == "time"
    k_ = find(C.t >= ONSET_TIME, 1);  if ~isempty(k_), k_on = k_; end
else
    for k = 1:height(C) - ONSET_RUN + 1
        if all(hit(k:k + ONSET_RUN - 1)), k_on = k; break; end
    end
end
after = false(height(C), 1);  if isfinite(k_on), after(k_on:end) = true; end
z2 = nan(height(C), 1);  R2 = nan(height(C), 1);
z2(after) = mod(C.psi(after) + PSI_OFFSET, 360);  R2(after) = PSI_SIGMA;
[C.track_psi, C.track_psi_sd, C.gated_th, C.gated_psi] = track_dir(C.t, C.theta, Rm, Q_RATE, GATE, GATE_RESET, z2, R2);
x_on = NaN;  if isfinite(k_on), x_on = x(k_on); end

%% CSV
werr = @(a, b) mod(a - b + 180, 360) - 180;          % signed error, wrapped to [-180, 180)
tcsv = C.t;  tcsv.TimeZone = '';  tcsv.Format = 'yyyy-MM-dd HH:mm:ss';
OUT = table(tcsv, round(C.Hm0, 2), round(C.Tm01, 1), round(C.theta, 1), round(C.track_psi, 1), ...
    round(werr(C.theta, C.ref_lp8), 1), round(werr(C.track_psi, C.ref_lp8), 1), ...
    round(C.U10, 1), round(C.wind_from, 0), round(C.ice, 2), round(C.coh, 3), ...
    'VariableNames', {'time_utc', 'Hm0_m', 'Tm01_s', 'dir_buoy_deg', 'dir_tracked_deg', ...
    'err_buoy_deg', 'err_tracked_deg', 'U10_mps', 'wind_from_deg', 'ice_conc', 'coherence'});
if ~exist(FIG.outdir, 'dir'), mkdir(FIG.outdir); end
csvfile = fullfile(FIG.outdir, [tag '_results.csv']);
writetable(OUT, csvfile);
fprintf('CSV: %s (%d rows)\n', csvfile, height(OUT));

%% F1: tracked direction 
fig = mkfig([1100 520]);
ax = axes(fig);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
wl(ax, x, C.ref_lp8, '--', [0 0 0], 'LP8 Mean Direction', FIG);
wl(ax, x, C.psi, '-', cg, '\psi (GNSS Heading Proxy)', FIG);
wl(ax, x, C.track, '-', cb, '\theta_{alg} (Earth Frame)', FIG);
wl(ax, x, C.track_psi, '-', c2, '\theta_{final} (With Heading)', FIG);
gt = C.gated_th == 1;  gp = C.gated_psi == 1;                 % gate counts (console only)
dirax(ax, x);
legend(ax, 'Location', 'northwest');
savepng(fig, FIG, FIG.outdir, [tag '_track_proxy']);

%% F2: onset rule
fig = mkfig([1100 600]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
sh = true(height(C), 1);                                       % deployment/recovery transients out
sh(1:min(F2_DROP_START, end)) = false;  sh(max(end - F2_DROP_END + 1, 1):end) = false;
a1 = nexttile(tl);  hold(a1, 'on');  grid(a1, 'on');  box(a1, 'on');
plot(a1, x(sh), C.yaw_max(sh), '-', 'Color', cb, 'LineWidth', FIG.lw);
yline(a1, YAW_ONSET, ':', sprintf('Threshold %g deg/s', YAW_ONSET), 'Color', cg, ...
    'LabelHorizontalAlignment', 'left', 'FontSize', FIG.fs - 2);
ylim(a1, [0 Inf]);
ylabel(a1, 'Max |Yaw Rate| (deg/s)');
title(a1, sprintf('%s: (a) Yaw Rate, (b) Kurtosis', BUOY_ID), 'FontWeight', 'normal');
a2 = nexttile(tl);  hold(a2, 'on');  grid(a2, 'on');  box(a2, 'on');
plot(a2, x(sh), C.kurt(sh), '-', 'Color', cb, 'LineWidth', FIG.lw);
set(a2, 'YScale', 'log');
yline(a2, KURT_GAUSS, '-', 'Gaussian (3)', 'Color', cg, 'LabelHorizontalAlignment', 'left', 'FontSize', FIG.fs - 2);
yline(a2, KURT_ONSET, ':', sprintf('Threshold %g', KURT_ONSET), 'Color', cg, ...
    'LabelHorizontalAlignment', 'right', 'FontSize', FIG.fs - 2);
ylabel(a2, 'Kurtosis of a_z > 2 Hz (-)');
xfin([a1 a2], x(sh));
savepng(fig, FIG, FIG.outdir, [tag '_onset_rule']);

%% F5: heave slope coherence and direction error 
% (a) coherence over the swell region (body frame, before tracking) (b) |error| of the
% per-packet direction against LP8. Same packets as F2 
fig = mkfig([1100 600]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
a1 = nexttile(tl);  hold(a1, 'on');  grid(a1, 'on');  box(a1, 'on');
plot(a1, x(sh), C.coh(sh), '-', 'Color', cb, 'LineWidth', FIG.lw);
yline(a1, R1_VALID^2, ':', sprintf('Validity Limit %.2f (r_1 = %.2f)', R1_VALID^2, R1_VALID), 'Color', cg, ...
    'LabelHorizontalAlignment', 'left', 'FontSize', FIG.fs - 2);
ylim(a1, [0 1]);  ylabel(a1, 'Swell-Band Coherence \gamma^2 (-)');
title(a1, sprintf('%s: (a) Heave-Slope Coherence, (b) Direction Error', BUOY_ID), 'FontWeight', 'normal');
a2 = nexttile(tl);  hold(a2, 'on');  grid(a2, 'on');  box(a2, 'on');
plot(a2, x(sh), abs(werr(C.theta(sh), C.ref_lp8(sh))), '-', 'Color', cb, 'LineWidth', FIG.lw);
ylim(a2, [0 180]);  yticks(a2, 0:45:180);  ylabel(a2, 'Error |\theta_{alg} - LP8| (deg)');
if isfinite(x_on)
    xline(a1, x_on, '--', 'Color', cg);
    xline(a2, x_on, '--', 'Onset', 'Color', cg, 'FontSize', FIG.fs - 2);
end
xfin([a1 a2], x(sh));
savepng(fig, FIG, FIG.outdir, [tag '_coherence']);

%%  F3 and F4: per-packet spectrum figure and wind figure 
if PACKET_FIGS == "all"
    pick = find(~cellfun(@isempty, SP))';
elseif isempty(PACKET_FIG_TIMES)
    pick = [];
    b_ = find(v1 & ~after);  if ~isempty(b_), [~, j_] = max(C.r1(b_));  pick(end+1) = b_(j_); end
    a_ = find(after);        if ~isempty(a_), [~, j_] = max(C.Hm0(a_)); pick(end+1) = a_(j_); end
else
    pick = arrayfun(@(tq) nearest_idx(C.t, tq), PACKET_FIG_TIMES);
end
dirP = FIG.outdir;  dirW = FIG.outdir;
if PACKET_FIGS == "all", dirP = fullfile(FIG.outdir, [tag '_packets']);  dirW = fullfile(FIG.outdir, [tag '_wind']); end
Umax = 5 * ceil(max([C.U10; 5], [], 'omitnan') / 5);
reset_plot_defaults();          
for i = unique(pick(:)')
    S = SP{i};  if isempty(S), continue; end
    stamp = char(datetime(C.t(i), 'Format', 'ddHHmm'));
    vstr  = iff(S.valid == 1, "Valid", "Invalid");
    head  = sprintf('%s  %s UTC  (%s Packet)', BUOY_ID, char(datetime(C.t(i), 'Format', 'dd-MMM-yyyy HH:mm')), vstr);
    % F3: polar E(f,theta) and Czz (layout of the batchCheck packet figure) ----
    fig = mkfig([1100 480]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
    title(tl, {head, sprintf(['H_{m0} = %.2f m,   \\theta_{alg} = %.0f^\\circ,   \\theta_{final} = %.0f^\\circ,   ' ...
        'LP8 = %.0f^\\circ,   r_1 = %.2f'], S.Hm0, S.theta, C.track_psi(i), C.ref_lp8(i), S.r1)});
    ax = nexttile(tl);  hold(ax, 'on');
    E = max(S.Czz(:) - S.Nzz(:), 0) .* S.D;
    dirs = struct('deg', {S.theta, C.track_psi(i), C.ref_lp8(i)}, 'col', {cb, c2, [0 0 0]}, ...
        'ls', {'-', '-', '--'}, 'name', {'\theta_{alg} (SHARC)', '\theta_{final} (Tracked)', 'LP8'});
    polar_E(ax, S.freq(:), S.tg(:), E, POLAR_FMAX, dirs, cmap);
    ax = nexttile(tl);
    k = S.freq > 0;
    semilogy(ax, S.freq(k), S.Czz(k), '-', 'Color', cb, 'LineWidth', FIG.lw, 'DisplayName', 'Heave Spectrum C_{zz}');
    hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
    semilogy(ax, S.freq(k), S.Nzz(k), 'k--', 'LineWidth', FIG.lw, 'DisplayName', 'IMU Noise N_{zz}');
    if any(S.band)
        yl = ylim(ax);  fb = S.freq(S.band);  fb = [fb(1) fb(end)];
        patch(ax, [fb fliplr(fb)], [yl(1) yl(1) yl(2) yl(2)], cband, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
            'DisplayName', 'SNR Band');
        ylim(ax, yl);
    end
    xlim(ax, [0 0.5]);  xlabel(ax, 'Frequency (Hz)');  ylabel(ax, 'Spectral Density (m^2/Hz)');
    legend(ax, 'Location', 'northeast');
    savepng(fig, FIG, dirP, sprintf('%s_packet_%s', tag, stamp));
    close(fig);
    %  F4: wind speed over the deployment (green line = this packet) and wind direction 
    fig = mkfig([1100 480]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
    title(tl, {head, sprintf('Wind Speed U_{10} = %.1f m/s,   From %.0f^\\circ', C.U10(i), C.wind_from(i))});
    ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
    plot(ax, x, C.U10, '-', 'Color', [0 0 0], 'LineWidth', FIG.lw);
    xline(ax, x(i), '-', 'Color', GREEN, 'LineWidth', 1.5);
    plot(ax, x(i), C.U10(i), 'o', 'Color', GREEN, 'MarkerFaceColor', GREEN, 'MarkerSize', 6);
    ylim(ax, [0 Umax]);  xlim(ax, [min(x) - 0.05, max(x) + 0.05]);
    xticks(ax, floor(min(x)):max(1, ceil((max(x) - min(x))/8)):ceil(max(x)));
    xlabel(ax, 'Day of July 2022 (UTC)');  ylabel(ax, 'Wind Speed U_{10} (m/s)');
    ax = nexttile(tl);  hold(ax, 'on');
    wind_compass(ax, C.wind_from(i), C.U10(i), Umax, GREEN);
    savepng(fig, FIG, dirW, sprintf('%s_wind_%s', tag, stamp));
    close(fig);
end
fprintf('F3/F4: %d packets -> %s and %s\n', numel(unique(pick)), dirP, dirW);

%% CONSOLE DEBUG
s_on = "not found";  if isfinite(k_on), s_on = string(C.t(k_on)); end   % iff evaluates both branches
fprintf('\nonset (%s): %s\n', ONSET_MODE, s_on);
fprintf('after onset: error vs LP8 (circular mean, std, n)\n');
est = {'per-packet, no tracking', C.theta; 'tracked, SHARC body', C.track; ...
       'tracked with heading', C.track_psi; 'psi alone', mod(C.psi + PSI_OFFSET, 360)};
for e = 1:size(est, 1)
    [m1, s1, n1] = cstat(est{e,2}, C.ref_lp8, after);
    fprintf('  %-26s %7.1f %7.1f %5d\n', est{e,1}, m1, s1, n1);
end
fprintf('  psi rejected %d of %d, theta rejected %d\n', nnz(gp & after), nnz(after & isfinite(z2)), nnz(gt & after));
fprintf('valid %d of %d packets. Figures in %s/ (%s_*)\n', nnz(v1), height(C), FIG.outdir, tag);


%% LOCAL FUNCTIONS
function reset_plot_defaults() %changes frm original plots
    p = {'defaultAxesFontName','defaultAxesFontSize','defaultAxesLineWidth','defaultAxesBox', ...
         'defaultAxesXGrid','defaultAxesYGrid','defaultAxesGridAlpha','defaultAxesGridLineStyle', ...
         'defaultAxesTickDir','defaultAxesTickDirMode','defaultAxesLayer','defaultAxesColorOrder', ...
         'defaultLineLineWidth','defaultLineMarkerSize','defaultTextFontName','defaultTextFontSize', ...
         'defaultLegendFontName','defaultLegendFontSize','defaultLegendBox','defaultColorbarFontName', ...
         'defaultColorbarFontSize','defaultPolarAxesFontName','defaultPolarAxesFontSize'};
    for k = 1:numel(p)
        try, set(groot, p{k}, 'remove'); catch, end
    end
    set(groot, 'defaultFigureColor', 'w');
end

function set_plot_defaults(fs)
% Titles and axis labels use the same size as the ticks (multipliers 1), legends and
% colorbars 2 pt smaller.
    set(groot, 'defaultAxesFontSize', fs, 'defaultTextFontSize', fs, ...
        'defaultAxesLabelFontSizeMultiplier', 1, 'defaultAxesTitleFontSizeMultiplier', 1, ...
        'defaultAxesTitleFontWeight', 'normal', ...
        'defaultLegendFontSize', fs - 2, 'defaultColorbarFontSize', fs - 2);
end

function fig = mkfig(px)
% Figure of px = [width height] pixels, white, 
    fig = figure('Color', 'w', 'Position', [60 60 px], 'Visible', 'off');
end

function savepng(fig, FIG, folder, name)
% PNG only, FIG.dpi.
    if ~exist(folder, 'dir'), mkdir(folder); end
    exportgraphics(fig, fullfile(folder, [name '.png']), 'Resolution', FIG.dpi, 'BackgroundColor', 'white');
end

function x = tday(t)
% Day of July 2022 (UTC) as a number: 21.5 = 21 July 12:00.
    t.TimeZone = '';  x = days(t - datetime(2022, 7, 1)) + 1;
end

function dirax(ax, x)
    ylim(ax, [0 360]);  yticks(ax, 0:90:360);  ylabel(ax, 'Mean Wave Direction (deg)');
    timeax(ax, x, true);
end

function xfin(axs, x)
    linkaxes(axs, 'x');
    for k = 1:numel(axs), timeax(axs(k), x, k == numel(axs)); end
end

function timeax(ax, x, lab)
% Time axis in day of July: labelled ticks every 6 h (HH:mm, date under 00:00) and a
% vertical grid line every hour (minor grid). lab = false: no tick labels, no xlabel.
% The label step widens (6 h, 12 h, 1 d, 2 d, 5 d) so that at most ~10 labels are shown
% at the report font size; only the tick labels change, not the data.
    x = x(isfinite(x));
    if isempty(x), warning('timeax: no finite times to plot (no packets?).'); return; end
    x0 = floor(min(x)*24)/24;  x1 = ceil(max(x)*24)/24;
    if x1 <= x0, x1 = x0 + 1/24; end                 % one packet exactly on the hour
    xlim(ax, [x0 x1]);
    stp = [0.25 0.5 1 2 5];  st = stp(find((x1 - x0)./stp <= 10, 1));  if isempty(st), st = 10; end
    mj = ceil(x0/st)*st : st : x1;
    xticks(ax, mj);
    ax.XAxis.MinorTickValues = x0 : 1/24 : x1;
    ax.XMinorTick = 'on';  ax.XMinorGrid = 'on';
    ax.MinorGridLineStyle = '-';  ax.MinorGridAlpha = 0.06;
    if ~lab, ax.XTickLabel = []; return; end
    L = strings(size(mj));
    for k = 1:numel(mj)
        hh = mod(round(mod(mj(k), 1)*24), 24);
        if hh == 0 && st >= 1
            L(k) = sprintf('%d Jul', round(mj(k)));
        elseif hh == 0
            L(k) = sprintf('00:00\\newline%d Jul', round(mj(k)));
        else
            L(k) = sprintf('%02d:00', hh);
        end
    end
    ax.TickLabelInterpreter = 'tex';
    xticklabels(ax, L);
    xlabel(ax, 'Time (UTC), July 2022');
end

function wl(ax, t, y, spec, c, name, FIG)
% Direction line with a gap where it wraps through 0/360 (no vertical streaks).
    t = t(:);  y = y(:);  k = isfinite(y) & isfinite(t);  t = t(k);  y = y(k);
    if isempty(y), return; end
    pos = cumsum([1; 1 + (abs(diff(y)) > 180)]);
    Y = nan(pos(end), 1);  Y(pos) = y;  T = nan(pos(end), 1);  T(pos) = t;
    h = plot(ax, T, Y, spec, 'Color', c, 'LineWidth', FIG.lw, 'MarkerSize', FIG.ms, 'MarkerFaceColor', 'w');
    if strlength(name) > 0, h.DisplayName = name; else, h.HandleVisibility = 'off'; end
end

function bnd(ax, t, y, s, c)
% +/- 1 sigma band around a direction series (split where it wraps).
    t = t(:);  y = y(:);  s = s(:);  k = isfinite(y) & isfinite(s) & isfinite(t);
    t = t(k);  y = y(k);  s = s(k);
    if numel(y) < 2, return; end
    br = [0; find(abs(diff(y)) > 180); numel(y)];
    for q = 1:numel(br) - 1
        i = br(q)+1:br(q+1);
        if numel(i) < 2, continue; end
        fill(ax, [t(i); flipud(t(i))], [y(i) - s(i); flipud(y(i) + s(i))], c, 'FaceAlpha', 0.15, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

function clim_(ax, lim)
    if exist('clim', 'file'), clim(ax, lim); else, caxis(ax, lim); end %#ok<CAXIS>
end

function k = nearest_idx(t, tq)
    t.TimeZone = '';  tq.TimeZone = '';  [~, k] = min(abs(t - tq));
end

function v = iff(c, a, b)
    if c, v = a; else, v = b; end
end

function [mu, sd, n] = cstat(est, ref, use)
% Circular mean and std (deg) of est - ref over the packets in use.
    e = mod(est - ref + 180, 360) - 180;  e = e(use & isfinite(e));  n = numel(e);
    mu = NaN;  sd = NaN;
    if n == 0, return; end
    z = mean(exp(1i*deg2rad(e)));
    mu = rad2deg(angle(z));  sd = rad2deg(sqrt(-2*log(min(max(abs(z), eps), 1 - eps))));
end

function polar_E(ax, f, th, E, fmax, dirs, cmap)
% Directional spectrum as filled contours on a compass rose
%   f (Hz), th direction FROM (deg, clockwise from N), E(f, th); dirs: .deg .col .ls .name
    k  = f > 0 & f <= fmax;  f = f(k);  E = E(k, :);
    th = [th(:); th(1) + 360];  E = E(:, [1:end 1]);
    E  = E / max(E(:), [], 'omitnan');
    [T, F] = meshgrid(th, f);  X = F .* sind(T);  Y = F .* cosd(T);
    axis(ax, 'equal');  axis(ax, 'off');
    gc = [0.80 0.80 0.80];  rings = 0.05:0.05:fmax;
    for rr = rings
        plot(ax, rr*sind(0:2:360), rr*cosd(0:2:360), '-', 'Color', gc, 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
    for a = 0:30:330
        plot(ax, [0 fmax*sind(a)], [0 fmax*cosd(a)], '-', 'Color', gc, 'LineWidth', 0.5, 'HandleVisibility', 'off');
        text(ax, 1.10*fmax*sind(a), 1.10*fmax*cosd(a), sprintf('%d', a), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    for rr = rings
        text(ax, rr*sind(15), rr*cosd(15), sprintf('%.2f', rr), 'FontSize', 8, ...
            'Color', [0.35 0.35 0.35], 'VerticalAlignment', 'bottom');
    end
    contourf(ax, X, Y, E, 0.1:0.1:1, 'LineStyle', 'none', 'HandleVisibility', 'off');
    contour(ax, X, Y, E, [0.02 0.05 0.1], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    colormap(ax, cmap);  clim_(ax, [0 1]);
    cbh = colorbar(ax);  cbh.Label.String = 'Normalised Energy E(f,\theta) (-)';
    for q = 1:numel(dirs)
        if isfinite(dirs(q).deg)
            plot(ax, [0 fmax*sind(dirs(q).deg)], [0 fmax*cosd(dirs(q).deg)], dirs(q).ls, ...
                'Color', dirs(q).col, 'LineWidth', 1.5, 'DisplayName', ...
                sprintf('%s %.0f^\\circ', dirs(q).name, dirs(q).deg));
        end
    end
    if any(arrayfun(@(d) isfinite(d.deg), dirs)), legend(ax, 'Location', 'southoutside'); end
    xlim(ax, 1.18*fmax*[-1 1]);  ylim(ax, 1.18*fmax*[-1 1]);
end

function wind_compass(ax, wdir, U, rmax, col)
% One packet's wind on a compass: line from the centre towards the direction the wind comes
% FROM (clockwise from N), length = U10 on rings every 5 m/s; N, E, S, W marked.
    axis(ax, 'equal');  axis(ax, 'off');  gc = [0.80 0.80 0.80];
    for rr = 5:5:rmax
        plot(ax, rr*sind(0:2:360), rr*cosd(0:2:360), '-', 'Color', gc, 'LineWidth', 0.5);
        text(ax, rr*sind(20), rr*cosd(20), sprintf('%d m/s', rr), 'FontSize', 8, 'Color', [0.35 0.35 0.35], ...
            'VerticalAlignment', 'bottom');
    end
    for a = 0:30:330
        plot(ax, [0 rmax*sind(a)], [0 rmax*cosd(a)], '-', 'Color', gc, 'LineWidth', 0.5);
    end
    lab = {'N', 'E', 'S', 'W'};
    for q = 1:4
        a = (q-1)*90;
        text(ax, 1.10*rmax*sind(a), 1.10*rmax*cosd(a), lab{q}, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontWeight', 'bold');
    end
    if isfinite(wdir) && isfinite(U)
        plot(ax, [0 U*sind(wdir)], [0 U*cosd(wdir)], '-', 'Color', col, 'LineWidth', 1.5);
        plot(ax, U*sind(wdir), U*cosd(wdir), 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 5);
    else
        text(ax, 0, 0, 'No ERA5 Wind', 'HorizontalAlignment', 'center');
    end
    xlim(ax, 1.18*rmax*[-1 1]);  ylim(ax, 1.18*rmax*[-1 1]);
end

function [xs, sd, gated, gated2] = track_dir(t, z, R, q, gate, nreset, z2, R2)
% Circular random-walk Kalman filter + RTS smoother (degrees). z: measurements (NaN = none),
% R: measurement std, q: process std per sqrt(hour). Optional z2, R2: a second measurement
% (heading proxy), gated the same way, applied after z. Returns smoothed x, std, gate flags.
    n = numel(z);  xs = nan(n,1);  sd = nan(n,1);  gated = zeros(n,1);  gated2 = zeros(n,1);
    if nargin < 7, z2 = nan(n,1);  R2 = nan(n,1); end
    w = @(a) mod(a + 180, 360) - 180;
    xf = nan(n,1);  Pf = nan(n,1);  xp = nan(n,1);  Pp = nan(n,1);
    k0 = find(isfinite(z) & R < max(R), 1);  if isempty(k0), k0 = find(isfinite(z), 1); end
    if isempty(k0), return; end
    x = z(k0);  P = R(k0)^2;  rej = [];
    for k = k0:n
        if k > k0, P = P + q^2 * max(hours(t(k) - t(k-1)), 0); end
        xp(k) = x;  Pp(k) = P;
        if isfinite(z(k))
            nu = w(z(k) - x);  Sk = P + R(k)^2;
            if abs(nu) <= gate*sqrt(Sk) || k == k0
                K = P / Sk;  x = x + K*nu;  P = (1 - K)*P;  rej = [];
            else
                gated(k) = 1;  rej(end+1) = k;                                  %#ok<AGROW>
                if numel(rej) >= nreset
                    zr = z(rej(end-nreset+1:end));
                    c = mean(exp(1i*deg2rad(zr)));
                    if rad2deg(sqrt(-2*log(max(abs(c), eps)))) < 30
                        x = rad2deg(angle(c));  P = mean(R(rej(end-nreset+1:end)).^2);
                        gated(rej(end-nreset+1:end)) = 0;  rej = [];
                    end
                end
            end
        end
        if isfinite(z2(k)) && isfinite(R2(k))
            nu = w(z2(k) - x);  Sk = P + R2(k)^2;
            if abs(nu) <= gate*sqrt(Sk)
                K = P / Sk;  x = x + K*nu;  P = (1 - K)*P;
            else
                gated2(k) = 1;
            end
        end
        xf(k) = x;  Pf(k) = P;
    end
    xs(n) = xf(n);  Ps = Pf(n);  sd(n) = sqrt(Ps);
    for k = n-1:-1:k0
        G = Pf(k) / Pp(k+1);
        xs(k) = xf(k) + G * w(xs(k+1) - xp(k+1));
        Ps = Pf(k) + G^2 * (Ps - Pp(k+1));
        sd(k) = sqrt(max(Ps, 0));
    end
    xs = mod(xs, 360);
end

function m = lp8_ref(R, tmid, win_min)
% LP8 Mdir: time weighted circular mean of the reference hours within +/- win_min.
    m = NaN;
    if isempty(R), return; end
    tmid.TimeZone = '';
    dtm = abs(minutes(R.TimeUTC - tmid));  sel = find(dtm <= win_min);
    a = [];  w = [];
    for j = sel'
        r = R(j,:);  mj = [];
        if ismember('OCR_text', r.Properties.VariableNames) && ~ismissing(r.OCR_text) && strlength(r.OCR_text) > 0
            seg = split(r.OCR_text, "|");  seg = seg(contains(seg, "Mdir"));
            if ~isempty(seg)
                rhs = extractAfter(seg(1), "=");
                if ~ismissing(rhs), mj = mod(str2double(regexp(rhs, '-?\d{1,3}(?:\.\d+)?', 'match')), 360); end
            end
        end
        if isempty(mj), mj = [r.Mdir_1, r.Mdir_2];  mj = mj(~isnan(mj)); end
        a = [a, mj];  w = [w, (1 - dtm(j)/win_min + eps)*ones(size(mj))];   %#ok<AGROW>
    end
    if ~isempty(a), m = mod(atan2d(sum(w.*sind(a)), sum(w.*cosd(a))), 360); end
end

function [la, lo] = pos_at(G, tm)
    la = NaN;  lo = NaN;
    if isempty(G), return; end
    t = G.DateTimeUTC;
    if ~isdatetime(t), t = datetime(string(t), 'InputFormat','dd-MMM-yyyy HH:mm:ss'); end
    t.TimeZone = '';  tm.TimeZone = '';
    ok = isfinite(G.Latitude) & isfinite(G.Longitude) & G.Latitude ~= 0;
    t = t(ok);  [t, iu] = unique(t);
    if numel(t) < 2, return; end
    a = G.Latitude(ok);  b = G.Longitude(ok);
    la = interp1(t, a(iu), tm, 'nearest', NaN);  lo = interp1(t, b(iu), tm, 'nearest', NaN);
end

function E = era5_load(fn, vars)
% lat, lon, time and the listed variables as lon x lat x time ([] if the file is absent).
    E = [];
    if ~isfile(fn), warning('%s not found: no wind in F4.', fn); return; end
    try
        info = ncinfo(fn);  vn = {info.Variables.Name};
        E.la = double(ncread(fn, vn{find(ismember(vn, {'latitude','lat'}), 1)}));
        E.lo = double(ncread(fn, vn{find(ismember(vn, {'longitude','lon'}), 1)}));
        tn = vn{find(ismember(vn, {'valid_time','time'}), 1)};
        tv = double(ncread(fn, tn));  u = ncreadatt(fn, tn, 'units');
        tok = regexp(u, '^\s*(\w+)\s+since\s+(\d{4}-\d{1,2}-\d{1,2})', 'tokens', 'once');
        t0 = datetime(tok{2}, 'InputFormat', 'yyyy-M-d');
        switch lower(tok{1})
            case {'seconds','second','s'}, E.t = t0 + seconds(tv);
            case {'hours','hour','h'},     E.t = t0 + hours(tv);
            case {'days','day'},           E.t = t0 + days(tv);
            otherwise,                     E.t = t0 + minutes(tv);
        end
        for k = 1:numel(vars)
            vi = info.Variables(strcmp(vn, vars{k}));
            if isempty(vi), warning('%s: %s not found.', fn, vars{k}); continue; end
            E.(vars{k}) = llt(double(ncread(fn, vars{k})), {vi.Dimensions.Name});
        end
    catch ME
        warning('%s not loaded: %s', fn, ME.message);  E = [];
    end
end

function X = llt(X, dims)
% Permute to lon x lat x time; extra dimensions merged (first finite value).
    il = find(ismember(dims, {'longitude','lon'}));  ia = find(ismember(dims, {'latitude','lat'}));
    it = find(ismember(dims, {'valid_time','time'}));  other = setdiff(1:numel(dims), [il ia it]);
    if isempty(il) || isempty(ia) || isempty(it), X = squeeze(X); return; end
    X = permute(X, [il ia it other]);
    if ~isempty(other)
        sz = size(X);  sz(end+1:3) = 1;  X = reshape(X, sz(1), sz(2), sz(3), []);
        Y = X(:,:,:,1);
        for k = 2:size(X, 4), m = isnan(Y);  Xk = X(:,:,:,k);  Y(m) = Xk(m); end
        X = Y;
    end
end

function v = era5_at(E, name, t, la, lo)
% Nearest cell, nearest hour (<= 1 h); NaN otherwise.
    v = NaN;
    if isempty(E) || ~isfield(E, name) || ~isfinite(la) || ~isfinite(lo) || isnat(t), return; end
    t.TimeZone = '';
    x = lo;  if max(E.lo) > 180, x = mod(x, 360); else, x = mod(x + 180, 360) - 180; end
    [~, io] = min(abs(E.lo - x));  [~, ia] = min(abs(E.la - la));  [dt, it] = min(abs(E.t - t));
    if hours(dt) <= 1, v = E.(name)(io, ia, it); end
end