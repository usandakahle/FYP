%% batch_run_sharc_packets.m
%
% 1) PIPELINE: runs run_sharc_dsp on every IMU burst in a window. Only SHARC data and
%    SHARC's own GNSS yaw go into the chain -- no reference data.
% 2) INTERNAL QC (no references): gyro vs GNSS heading consistency and the
%    burst-to-burst step of the earth-frame direction.
% 3) VALIDATION (batch only): matches each burst to the LP8 figure directions, the
%    SB04 Hm0/Tm01 reference CSV and LP4, and exports everything with the design
%    settings to Excel (one sheet per run + a summary row).
%
% Matching is by TIME, never by packet number.

clc;

%% ------------------------------------------------------------
% PREREQUISITE: Results.YawEstimate_deg (SHARC GNSS yaw)
% ------------------------------------------------------------
if ~exist('Results','var') || ~ismember('YawEstimate_deg', Results.Properties.VariableNames)
    if isfile('SB04_yaw_state.mat')
        fprintf('Loading SB04_yaw_state.mat checkpoint.\n');
        load('SB04_yaw_state.mat', 'Results');
    else
        error('Results.YawEstimate_deg missing and no SB04_yaw_state.mat checkpoint.');
    end
end

%% ------------------------------------------------------------
% CONFIG
% ------------------------------------------------------------
PACKET_DIR            = "SB04/sb04_packets";
HEAVE_IS_DISPLACEMENT = false;
SCALE_HEAVE           = 9.80665 / 8192;   % ICM20649 +/-4 g, 16-bit: confirm FS setting
SCALE_ANGLE_DEG       = 1;                % unused by the chain
SIGN_QUAD             = +1;
SCALE_GYRO_RATE       = 1/65.5;           % must match run_sharc_dsp

windowStart    = datetime(2022,7,21,15,0,0,  'TimeZone','UTC');
windowEnd      = datetime(2022,7,22,01,59,0, 'TimeZone','UTC');
maxTimeGap_sec = 20*60;

% ---- HEADING ----
HEADING_MODE  = "in_chain";  % "in_chain": gyro yaw anchored to the GNSS fix, slopes rotated per sample
                             % "posthoc" : body-frame theta + GNSS yaw at the burst midpoint
YAW_SIGN      = 1;           % +1 if the Yaw gyro count is clockwise-positive viewed from above
YAW_BIAS_MODE = "per_burst"; % "per_burst": remove each burst's mean yaw rate (no net rotation in a burst)
                             % "prior"    : median-of-bursts bias (only if the gyro bias is stable)

% ---- DESIGN SETTINGS (exported with every run) ----
ALPHA_C       = 0.999;   % complementary filter; fc = (1-alpha)/(2*pi*dt)
TAU           = 100;     % RC high-pass (s); fc = 1/(2*pi*tau)
LPF_CUTOFF_HZ = 1.12;
LPF_ORDER     = 64;
M_DEC         = 32;
NFFT          = 512;
F_MIN         = 0.03;    % lower search limit (auto) or fixed cut (fixed)
F_MIN_MODE    = "auto";  % "auto": per-burst spectral minimum between F_MIN and 0.8*fp
F_MAX_BAND    = 0.3;
RUN_NOTE      = "in-chain heading, auto f_min";

% ---- VALIDATION REFERENCES (never passed to the chain) ----
REF_XLSX        = "reference_directions.xlsx";  % LP8 Mdir/Pdir from extract_reference_directions.m
REF_DIR_WIN_MIN = 60;                           % LP8 hours within +/- this of the burst midpoint, time-weighted
REF_PDIR_MAXSPR = 25;                           % LP8 Pdir values spread above this (deg) = bimodal hour
REF_HS_CSV      = "SB04_Hm0_reference.csv";     % DateTimeUTC,Hm0,Tm01,Hm0_upper,Hm0_lower
REF_HS_MAX_DT_MIN = 5;                          % matched to burst START
LP4_FILE        = "LP4.csv";
LP4_PAD_MIN     = 11;
ERA5_NC         = "ERA5Ice.nc";   % needs siconc; u10/v10 optional (wind-sea check)
ERA5_MAX_DT_H   = 1;
SIC_CLASS_EDGES = [0.15 0.80];      % <15% open, 15-80% MIZ, >80% pack (WMO convention)

% ---- UNCERTAINTY BUDGET (1-sigma, deg) -- justify each in the report ----
SIGMA_YAW_GNSS  = 10;   % GNSS-derived yaw error: set from the yaw filter's residuals

OUT_XLSX        = "SB04_direction_tuning.xlsx";

opts = struct('make_plots', true, 'verbose', true, 'station', 'SB04', ...
    'tau', TAU, 'M', M_DEC, 'lpf_cutoff_hz', LPF_CUTOFF_HZ, 'lpf_order', LPF_ORDER, ...
    'Nfft', NFFT, 'f_min', F_MIN, 'f_min_mode', F_MIN_MODE, 'f_max_band', F_MAX_BAND, ...
    'yaw_sign', YAW_SIGN, 'gnss_max_gap_s', maxTimeGap_sec, ...
    'regime_aware', true, 'dir_estimator', 'mem', 'despike', true, ...
    'regime_thr', struct('n_ice', 7, 'n_open', 5, 'hf_open', 0.15, 'kurt_miz', 6, ...
                         'hyst', 1, 'persist', 2), ...                       % calibrate vs SIC
    'yaw_bias_mode', YAW_BIAS_MODE, 'plot_estimators', {{'dfs','mem'}}, ...
    'qc_kurt_max', 100, 'qc_yaw_rms_max', 0.5, 'fp_min', 0.045);

%% ------------------------------------------------------------
% GNSS YAW (SHARC's own) and FILES IN WINDOW
% ------------------------------------------------------------
G = table(Results.TimeUTC, Results.YawEstimate_deg, 'VariableNames', {'TimeUTC','YawEstimate_deg'});
if isempty(G.TimeUTC.TimeZone), G.TimeUTC.TimeZone = 'UTC'; end
G = sortrows(G(~isnan(G.YawEstimate_deg) & ~isnat(G.TimeUTC), :), 'TimeUTC');
tG = posixtime(G.TimeUTC);

files = dir(fullfile(PACKET_DIR, "packet_*.csv"));
if isempty(files), error('No packet files found in %s.', PACKET_DIR); end
file_start = NaT(numel(files),1, 'TimeZone','UTC');
for k = 1:numel(files)
    tok = regexp(files(k).name, 'packet_(\d+)_(.+)\.csv', 'tokens', 'once');
    if ~isempty(tok)
        file_start(k) = datetime(tok{2}, 'InputFormat','yyyy-MM-dd_HH-mm-ss', 'TimeZone','UTC');
    end
end
selected_files = files(file_start >= windowStart & file_start <= windowEnd);
nP = numel(selected_files);
fprintf('Found %d packet files in [%s, %s].\n', nP, windowStart, windowEnd);
if nP == 0, error('No files in the requested window.'); end

%% Yaw-rate bias prior: median of per-burst mean gyro yaw rate (assumes no systematic spin)
pre_rate = nan(nP,1);
for f = 1:nP
    op = detectImportOptions(fullfile(selected_files(f).folder, selected_files(f).name));
    op.SelectedVariableNames = {'Yaw'};
    Ty = readtable(fullfile(selected_files(f).folder, selected_files(f).name), op);
    pre_rate(f) = YAW_SIGN * SCALE_GYRO_RATE * mean(double(Ty.Yaw), 'omitnan');
end
yaw_bias_prior = median(pre_rate, 'omitnan');
fprintf('Yaw-rate bias prior: %+.4f deg/s (IQR %.4f)\n', yaw_bias_prior, iqr(pre_rate));

%% ------------------------------------------------------------
% PIPELINE: run the chain on each burst
% ------------------------------------------------------------
rows = struct([]);  SET = [];
reg_state = struct('regime', "", 'pending', "", 'count', 0);   % carried burst to burst (causal)
for f = 1:nP
    fpath = fullfile(selected_files(f).folder, selected_files(f).name);
    Tw = readtable(fpath);
    if ~isdatetime(Tw.DateTimeUTC)
        Tw.DateTimeUTC = datetime(Tw.DateTimeUTC, 'InputFormat','dd-MMM-yyyy HH:mm:ss', 'TimeZone','UTC');
    elseif isempty(Tw.DateTimeUTC.TimeZone)
        Tw.DateTimeUTC.TimeZone = 'UTC';
    end
    label = sprintf('%s (%d rows)', selected_files(f).name, height(Tw));
    if height(Tw) < 2, warning('%s: insufficient data, skipping.', label); continue; end
    burstStart = min(Tw.DateTimeUTC);  burstEnd = max(Tw.DateTimeUTC);
    burstMid   = burstStart + (burstEnd - burstStart)/2;
    fprintf('\n=== [%d/%d] %s ===\nBurst: %s to %s\n', f, nP, label, burstStart, burstEnd);

    [gap_s, iy] = min(abs(tG - posixtime(burstMid)));
    yaw_mid = NaN;  if gap_s <= maxTimeGap_sec, yaw_mid = G.YawEstimate_deg(iy); end

    opts_i = opts;
    opts_i.yaw_rate_bias_dps = yaw_bias_prior;
    opts_i.regime_state = reg_state;
    if HEADING_MODE == "in_chain"
        opts_i.gnss_yaw = G;            % chain builds psi(t) and rotates slopes per sample
    else
        opts_i.yaw_deg = yaw_mid;       % figures only: body + yaw
    end

    try
        [theta, Hm0, Tm01, ~, pd] = run_sharc_dsp(Tw, label, HEAVE_IS_DISPLACEMENT, SCALE_HEAVE, ...
            SCALE_ANGLE_DEG, SIGN_QUAD, ALPHA_C, opts_i);
    catch ME
        warning('%s failed: %s', selected_files(f).name, ME.message);
        continue;
    end
    SET = pd.settings;
    reg_state = pd.regime_state;

    r = struct();
    r.FileLabel   = string(selected_files(f).name);
    r.BurstStart  = burstStart;  r.BurstMid = burstMid;  r.BurstEnd = burstEnd;
    r.Hm0 = Hm0;  r.Tm01 = Tm01;  r.fp_Hz = pd.fp;
    r.Tm01_over_Tp = pd.Tm01_over_Tp;  r.f_lo_Hz = pd.f_lo_Hz;  r.m0_frac_below_0p05 = pd.m0_frac_below_0p05;
    r.R1 = pd.R1_summary;  r.slope_snr = pd.slope_snr;  r.spread_peak_deg = pd.spread_peak;
    r.heading_mode = string(pd.heading.mode);  r.n_fix = pd.heading.n_fix;
    r.yaw_mid = yaw_mid;  r.yaw_gap_min = gap_s/60;
    r.gyro_dyaw_deg = pd.gyro_dyaw_deg;
    r.qc_ok = pd.qc_ok;  r.regime_raw = string(pd.regime_raw);
    r.regime = string(pd.regime);  r.tail_n = pd.tail_n;  r.hf_frac = pd.hf_frac;
    r.kurt_hf = pd.kurt_hf;  r.n_spikes = pd.n_spikes;  r.yaw_rms_wave = pd.yaw_rms_wave;
    r.f_hi_Hz = pd.f_hi_Hz;  r.sigma_stat_deg = pd.theta_jack_sd;  r.Hm0_sd = pd.hm0_sd;
    r.n_peaks = pd.part.n_peaks;  r.e_wind = pd.part.e_wind;
    r.f_swell = pd.part.f_swell;  r.f_wind = pd.part.f_wind;
    yrot = 0;  if pd.part_frame == "body", yrot = yaw_mid; end    % partition dirs to earth frame
    r.Mdir_swell = mod(pd.part.mdir_swell + yrot, 360);
    r.Mdir_wind  = mod(pd.part.mdir_wind  + yrot, 360);
    r.Pdir_MEM   = mod(pd.pdir_mem + yrot, 360);
    r.theta_body = pd.theta_body;  r.theta_body_w = pd.theta_body_w;
    r.pdir_body  = pd.pdir_body;   r.pdir_body_w  = pd.pdir_body_w;
    if HEADING_MODE == "in_chain"
        r.Mdir_B = pd.theta_earth_stbd;  r.Mdir_W = pd.theta_earth_stbd_w;
        r.Pdir_B = pd.pdir_earth_stbd;   r.Pdir_W = pd.pdir_earth_stbd_w;
    else
        r.Mdir_B = mod(pd.theta_body  + yaw_mid, 360);  r.Mdir_W = mod(pd.theta_body_w + yaw_mid, 360);
        r.Pdir_B = mod(pd.pdir_body   + yaw_mid, 360);  r.Pdir_W = mod(pd.pdir_body_w  + yaw_mid, 360);
    end
    if isempty(rows), rows = r; else, rows(end+1) = r; end %#ok<AGROW>
end
if isempty(rows), error('No bursts processed successfully.'); end
DSP = struct2table(rows);
fprintf('\nProcessed %d of %d bursts.\n', height(DSP), nP);

%% ------------------------------------------------------------
% INTERNAL QC (no references)
% ------------------------------------------------------------
% (a) Burst-to-burst step of the earth-frame direction. Wave direction normally turns
%     slowly (~<15 deg per 40 min); large steps point at the heading, not the DSP.
DSP.dir_step_deg = [NaN; wrap180(diff(DSP.Mdir_B))];
DSP.body_step_deg = [NaN; wrap180(diff(DSP.theta_body))];
DSP.yaw_step_deg  = [NaN; wrap180(diff(DSP.yaw_mid))];
% (b) Rotation the gyro sees inside each burst vs the GNSS yaw change to the next fix.
%     If the gyro shows tens of degrees and GNSS yaw barely moves, GNSS yaw is not hull heading.
fprintf('\n--- INTERNAL HEADING QC ---\n');
fprintf('median |earth-frame direction step| = %.1f deg, max %.1f deg\n', ...
    median(abs(DSP.dir_step_deg), 'omitnan'), max(abs(DSP.dir_step_deg)));
fprintf('median |body-frame step| = %.1f deg vs median |GNSS yaw step| = %.1f deg\n', ...
    median(abs(DSP.body_step_deg), 'omitnan'), median(abs(DSP.yaw_step_deg), 'omitnan'));
fprintf('median |gyro rotation within burst| = %.1f deg (max %.1f)\n', ...
    median(abs(DSP.gyro_dyaw_deg), 'omitnan'), max(abs(DSP.gyro_dyaw_deg)));
ok = ~isnan(DSP.body_step_deg) & ~isnan(DSP.yaw_step_deg);
if nnz(ok) >= 4
    c = corrcoef(DSP.body_step_deg(ok), -DSP.yaw_step_deg(ok));
    fprintf(['corr(body step, -GNSS yaw step) = %.2f  (near +1: hull rotation is captured by ' ...
        'GNSS yaw; near 0: it is not)\n'], c(1,2));
end

%% ------------------------------------------------------------
% VALIDATION (batch only): ERA5 ice concentration, ice-edge normal, wind
% ------------------------------------------------------------
DSP = add_era5(DSP, Results, ERA5_NC, ERA5_MAX_DT_H, maxTimeGap_sec);
DSP.SIC_class = repmat("unknown", height(DSP), 1);
DSP.SIC_class(DSP.SIC <  SIC_CLASS_EDGES(1)) = "open";
DSP.SIC_class(DSP.SIC >= SIC_CLASS_EDGES(1) & DSP.SIC < SIC_CLASS_EDGES(2)) = "miz";
DSP.SIC_class(DSP.SIC >= SIC_CLASS_EDGES(2)) = "ice";

% Uncertainty budget (pipeline quantities + stated heading error)
sig_bias = iqr(pre_rate) / 1.349;                           % spread of per-burst mean yaw rate (deg/s)
T_half   = seconds(DSP.BurstEnd - DSP.BurstStart) / 2;
if HEADING_MODE == "in_chain" && YAW_BIAS_MODE == "per_burst"
    DSP.sigma_head_deg = SIGMA_YAW_GNSS * ones(height(DSP), 1);   % net in-burst rotation not resolvable
elseif HEADING_MODE == "in_chain"
    DSP.sigma_head_deg = sqrt(SIGMA_YAW_GNSS^2 + (sig_bias .* T_half).^2);
else                                                         % within-burst rotation not removed
    DSP.sigma_head_deg = sqrt(SIGMA_YAW_GNSS^2 + (DSP.gyro_dyaw_deg/2).^2);
end

%% ------------------------------------------------------------
% EXPORT
% ------------------------------------------------------------
export_tuning(DSP, SET, REF_XLSX, REF_DIR_WIN_MIN, REF_PDIR_MAXSPR, REF_HS_CSV, REF_HS_MAX_DT_MIN, ...
    LP4_FILE, LP4_PAD_MIN, OUT_XLSX, RUN_NOTE, HEADING_MODE, ...
    sprintf('%s to %s', windowStart, windowEnd));

%% ========================================================================
%% Local functions
%% ========================================================================
function export_tuning(DSP, SET, REF_XLSX, win_min, pdir_maxspr, REF_HS_CSV, hs_dt_min, ...
        LP4_FILE, lp4_pad, OUT_XLSX, note, heading_mode, window_txt)
    n = height(DSP);
    e180 = @(a, b) wrap180(a - b);

    E = table();
    E.BurstMidUTC    = DSP.BurstMid;
    E.File           = DSP.FileLabel;
    E.HeadingMode    = repmat(string(heading_mode), n, 1);
    E.HeadingFix     = DSP.heading_mode;
    E.Yaw_mid_deg    = DSP.yaw_mid;
    E.gyro_dyaw_deg  = DSP.gyro_dyaw_deg;
    E.dir_step_deg   = DSP.dir_step_deg;
    E.Hm0_m          = DSP.Hm0;
    E.Tm01_s         = DSP.Tm01;
    E.fp_Hz          = DSP.fp_Hz;
    E.Tm01_over_Tp   = DSP.Tm01_over_Tp;
    E.f_lo_Hz        = DSP.f_lo_Hz;
    E.m0_frac_below_0p05 = DSP.m0_frac_below_0p05;
    E.Mdir_SHARC_B   = DSP.Mdir_B;   E.Mdir_SHARC_W = DSP.Mdir_W;
    E.Pdir_SHARC_B   = DSP.Pdir_B;   E.Pdir_SHARC_W = DSP.Pdir_W;
    E.R1             = DSP.R1;
    E.spread_peak_deg = DSP.spread_peak_deg;
    E.Pdir_MEM       = DSP.Pdir_MEM;
    % regime (pipeline) and its features
    E.QC_ok          = DSP.qc_ok;
    E.Regime         = DSP.regime;
    E.Regime_raw     = DSP.regime_raw;
    E.theta_body_B   = DSP.theta_body;
    E.tail_n         = DSP.tail_n;
    E.hf_frac        = DSP.hf_frac;
    E.kurt_hf        = DSP.kurt_hf;
    E.n_spikes       = DSP.n_spikes;
    E.yaw_rms_wave   = DSP.yaw_rms_wave;
    E.f_hi_Hz        = DSP.f_hi_Hz;
    E.n_peaks        = DSP.n_peaks;
    E.f_swell_Hz     = DSP.f_swell;   E.f_wind_Hz = DSP.f_wind;   E.e_wind = DSP.e_wind;
    E.Mdir_swell     = DSP.Mdir_swell;   E.Mdir_wind = DSP.Mdir_wind;
    % ERA5 (validation only)
    E.SIC            = DSP.SIC;
    E.SIC_class      = DSP.SIC_class;
    E.ice_normal_from_deg = DSP.ice_normal_from_deg;
    E.wind_from_deg  = DSP.wind_from_deg;   E.wind_speed = DSP.wind_speed;
    E.err_wind_sea_vs_wind = wrap180(DSP.Mdir_wind - DSP.wind_from_deg);
    % uncertainty
    E.sigma_stat_deg = DSP.sigma_stat_deg;
    E.sigma_head_deg = DSP.sigma_head_deg;
    E.Hm0_sd_m       = DSP.Hm0_sd;

    % ---- LP8 directions (time-weighted over +/- win_min) ----
    E.Ref_Mdir = nan(n,1);  E.Ref_Pdir = nan(n,1);  E.Ref_Pdir_spread = nan(n,1);  E.Ref_hours = strings(n,1);
    if isfile(REF_XLSX)
        R = readtable(REF_XLSX, 'Sheet', 'reference', 'TextType', 'string');
        if ~isdatetime(R.TimeUTC), R.TimeUTC = datetime(R.TimeUTC); end
        R.TimeUTC.TimeZone = 'UTC';
        for i = 1:n
            dtm = abs(minutes(R.TimeUTC - DSP.BurstMid(i)));
            sel = find(dtm <= win_min);
            if isempty(sel), continue; end
            m = []; p = []; wm = []; wp = [];
            for j = sel'
                [mj, pj] = ref_values(R(j,:));
                w = 1 - dtm(j)/win_min + eps;
                m = [m, mj]; wm = [wm, w*ones(size(mj))]; %#ok<AGROW>
                p = [p, pj]; wp = [wp, w*ones(size(pj))]; %#ok<AGROW>
            end
            E.Ref_Mdir(i) = wcircmean(m, wm);
            E.Ref_Pdir(i) = wcircmean(p, wp);
            [~, E.Ref_Pdir_spread(i)] = circ_stats(wrap180(p - E.Ref_Pdir(i)));
            E.Ref_hours(i) = strjoin(string(datetime(R.TimeUTC(sel), 'Format', 'dd HH:mm')), ', ');
        end
    else
        warning('%s not found: direction validation skipped.', REF_XLSX);
    end
    E.err_Mdir_B = e180(E.Mdir_SHARC_B, E.Ref_Mdir);
    E.err_Mdir_W = e180(E.Mdir_SHARC_W, E.Ref_Mdir);
    E.err_Pdir_B = e180(E.Pdir_SHARC_B, E.Ref_Pdir);
    E.err_Pdir_W = e180(E.Pdir_SHARC_W, E.Ref_Pdir);
    E.Ref_Pdir_bimodal = E.Ref_Pdir_spread > pdir_maxspr;

    % LP8 spread within the window = reference uncertainty; total predicted sigma
    E.sigma_ref_deg  = nan(n,1);
    if isfile(REF_XLSX)
        for i = 1:n
            dtm = abs(minutes(R.TimeUTC - DSP.BurstMid(i)));  sel = find(dtm <= win_min);
            m = [];  for j = sel', m = [m, ref_values(R(j,:))]; end %#ok<AGROW>
            if numel(m) >= 2, [~, E.sigma_ref_deg(i)] = circ_stats(wrap180(m - wcircmean(m, ones(size(m))))); end
        end
    end
    E.sigma_total_deg = sqrt(E.sigma_stat_deg.^2 + E.sigma_head_deg.^2 + E.sigma_ref_deg.^2);
    E.within_2sigma   = abs(E.err_Mdir_B) <= 2*E.sigma_total_deg;
    % directional filtering by ice: in-ice direction should lean toward the ice-edge normal
    E.LP8_minus_normal = wrap180(E.Ref_Mdir - E.ice_normal_from_deg);

    % ---- SB04 Hm0/Tm01 reference CSV (matched to burst start) ----
    E.RefCSV_Hm0 = nan(n,1);  E.RefCSV_Tm01 = nan(n,1);
    if isfile(REF_HS_CSV)
        H = readtable(REF_HS_CSV);
        if ~isdatetime(H.DateTimeUTC)
            H.DateTimeUTC = datetime(H.DateTimeUTC, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
        end
        H.DateTimeUTC.TimeZone = 'UTC';
        H = H(H.Hm0 > 0 & ~isnan(H.Tm01), :);
        for i = 1:n
            [dt, k] = min(abs(H.DateTimeUTC - DSP.BurstStart(i)));
            if ~isempty(k) && minutes(dt) <= hs_dt_min
                E.RefCSV_Hm0(i) = H.Hm0(k);  E.RefCSV_Tm01(i) = H.Tm01(k);
            end
        end
    else
        warning('%s not found: Hm0 reference skipped.', REF_HS_CSV);
    end
    E.Hm0_ratio_CSV  = E.RefCSV_Hm0  ./ E.Hm0_m;
    E.Tm01_ratio_CSV = E.RefCSV_Tm01 ./ E.Tm01_s;

    % ---- LP4 (in-situ): Hm0, Tm01, Meandir over the burst window ----
    E.LP4_Hm0 = nan(n,1);  E.LP4_Tm01 = nan(n,1);  E.LP4_Meandir = nan(n,1);
    if isfile(LP4_FILE)
        L = load_lp4(LP4_FILE);
        for i = 1:n
            s = L.TimeUTC >= DSP.BurstStart(i) - minutes(lp4_pad) & L.TimeUTC <= DSP.BurstEnd(i) + minutes(lp4_pad);
            if ~any(s), continue; end
            E.LP4_Hm0(i)  = mean(L.Hm0(s), 'omitnan');
            E.LP4_Tm01(i) = mean(L.Tm01(s), 'omitnan');
            E.LP4_Meandir(i) = wcircmean(L.Meandir(s)', ones(1, nnz(s)));
        end
    end
    E.Hm0_ratio_LP4    = E.LP4_Hm0 ./ E.Hm0_m;
    E.err_Mdir_B_LP4   = e180(E.Mdir_SHARC_B, E.LP4_Meandir);

    if ~isempty(SET), E = [E, repmat(struct2table(SET), n, 1)]; end

    runTime = datetime('now', 'TimeZone', 'UTC');
    sheet = "run_" + string(datetime(runTime, 'Format', 'yyyyMMdd_HHmmss'));
    writetable(strip_tz(E), OUT_XLSX, 'Sheet', sheet);

    % ---- summary row ----
    S = table(strip_tz(runTime), string(note), string(window_txt), sheet, string(heading_mode), ...
        'VariableNames', {'RunUTC','Note','Window','Sheet','HeadingMode'});
    if ~isempty(SET), S = [S, struct2table(SET)]; end
    S.N = n;
    S.N_QC_fail = nnz(~E.QC_ok);
    E0 = E;                                   % sheet keeps every burst; stats use QC-passed only
    for c = ["err_Mdir_B","err_Mdir_W","err_Pdir_B","err_Pdir_W","err_Mdir_B_LP4","dir_step_deg"]
        E.(c)(~E.QC_ok) = NaN;
    end
    S.N_dir_matched = nnz(~isnan(E.err_Mdir_B));
    S = add_err(S, "Mdir_B", E.err_Mdir_B);
    S = add_err(S, "Mdir_W", E.err_Mdir_W);
    S = add_err(S, "Pdir_B", E.err_Pdir_B);
    S = add_err(S, "Pdir_B_unimodal", E.err_Pdir_B(~E.Ref_Pdir_bimodal));
    S = add_err(S, "Mdir_B_vsLP4", E.err_Mdir_B_LP4);
    S.dir_step_medabs  = median(abs(E.dir_step_deg), 'omitnan');
    S.Hm0_ratio_CSV_med = median(E.Hm0_ratio_CSV, 'omitnan');
    S.Hm0_ratio_CSV_iqr = iqr(E.Hm0_ratio_CSV(~isnan(E.Hm0_ratio_CSV)));
    S.Tm01_ratio_CSV_med = median(E.Tm01_ratio_CSV, 'omitnan');
    S.Hm0_ratio_LP4_med = median(E.Hm0_ratio_LP4, 'omitnan');
    S.Tm01_over_Tp_med  = median(E.Tm01_over_Tp, 'omitnan');
    S.f_lo_med          = median(E.f_lo_Hz, 'omitnan');
    S.R1_median         = median(E.R1, 'omitnan');
    v = ~isnan(E.err_Mdir_B) & ~isnan(E.sigma_total_deg);
    S.coverage_2sigma   = mean(E.within_2sigma(v));             % ~0.95 if the budget is right
    S.sigma_total_med   = median(E.sigma_total_deg(v));
    for g = ["ice","miz","open"]
        k = E.Regime == g;
        S.("N_" + g) = nnz(k);
        e = E.err_Mdir_B(k & ~isnan(E.err_Mdir_B));
        [~, S.("Mdir_cstd_" + g)] = circ_stats(e);
        S.("spread_med_" + g) = median(E.spread_peak_deg(k), 'omitnan');
    end
    % detector vs ERA5 SIC
    known = E.SIC_class ~= "unknown";
    S.regime_agree_SIC = mean(E.Regime(known) == E.SIC_class(known));
    S.regime_switches  = nnz(E.Regime(2:end) ~= E.Regime(1:end-1) & E.QC_ok(2:end) & E.QC_ok(1:end-1));
    S.raw_switches     = nnz(E.Regime_raw(2:end) ~= E.Regime_raw(1:end-1));
    ok = ~isnan(E.LP8_minus_normal) & ~isnan(E.err_Mdir_B) & E.Regime ~= "open";
    if nnz(ok) >= 4
        c = corrcoef(-E.LP8_minus_normal(ok), E.err_Mdir_B(ok));  S.corr_err_vs_edge = c(1,2);
    else
        S.corr_err_vs_edge = NaN;
    end

    write_summary(S, OUT_XLSX, 'summary_v3');     % rows merged by column NAME (v2 was misaligned)
    E = E0; %#ok<NASGU>
    fprintf('\nRegime (pipeline) vs ERA5 SIC class:\n');
    disp(crosstab_table(E.Regime, E.SIC_class));
    fprintf(['2-sigma coverage of |Mdir error| = %.0f%% (target ~95%%); corr(error, edge pull) = %.2f\n'], ...
        100*S.coverage_2sigma, S.corr_err_vs_edge);
    fprintf(['\n--- VALIDATION (%s) ---\nMdir vs LP8: mean %+.1f, circ std %.1f, median|e| %.1f deg (N=%d)\n' ...
        'Hm0 ref-CSV/SHARC median %.2f, Tm01 ratio %.3f; Hm0 LP4/SHARC median %.2f; Tm01/Tp median %.2f\n' ...
        'Wrote sheet %s in %s\n'], heading_mode, S.Mdir_B_mean, S.Mdir_B_cstd, S.Mdir_B_medabs, ...
        S.N_dir_matched, S.Hm0_ratio_CSV_med, S.Tm01_ratio_CSV_med, S.Hm0_ratio_LP4_med, ...
        S.Tm01_over_Tp_med, sheet, OUT_XLSX);
end

function S = add_err(S, name, e)
    e = e(~isnan(e));
    [mu, sd] = circ_stats(e);
    S.(name + "_mean") = mu;  S.(name + "_cstd") = sd;  S.(name + "_medabs") = median(abs(e));
end

function [m, p] = ref_values(r)
% All values in the OCR line (one per LP8 record in that hour), signs kept.
    m = []; p = [];
    if ismember('OCR_text', r.Properties.VariableNames) && ~ismissing(r.OCR_text) && strlength(r.OCR_text) > 0
        seg = split(r.OCR_text, "|");
        m = nums_after_eq(seg(contains(seg, "Mdir")));
        p = nums_after_eq(seg(contains(seg, "Pdir")));
    end
    if isempty(m), m = [r.Mdir_1, r.Mdir_2]; m = m(~isnan(m)); end
    if isempty(p), p = [r.Pdir_1, r.Pdir_2]; p = p(~isnan(p)); end
end

function v = nums_after_eq(seg)
    v = [];
    if isempty(seg), return; end
    rhs = extractAfter(seg(1), "=");
    if ismissing(rhs), return; end
    v = mod(str2double(regexp(rhs, '-?\d{1,3}(?:\.\d+)?', 'match')), 360);
end

function m = wcircmean(a, w)
    k = ~isnan(a);  a = a(k);  w = w(k);
    if isempty(a), m = NaN; return; end
    m = mod(atan2d(sum(w.*sind(a)), sum(w.*cosd(a))), 360);
end

function [mu, sd] = circ_stats(e)
    if isempty(e), mu = NaN; sd = NaN; return; end
    z  = mean(exp(1i*deg2rad(e)));
    mu = rad2deg(angle(z));
    sd = rad2deg(sqrt(-2*log(min(max(abs(z), eps), 1 - eps))));
end

function y = wrap180(x)
    y = mod(x + 180, 360) - 180;
end

function T = strip_tz(T)
    if isdatetime(T), T.TimeZone = ''; return; end
    for v = string(T.Properties.VariableNames)
        if isdatetime(T.(v)), T.(v).TimeZone = ''; end
    end
end

function L = load_lp4(fn)
    op = detectImportOptions(fn, 'VariableNamingRule', 'preserve');
    T  = readtable(fn, op);
    vn = lower(strtrim(string(T.Properties.VariableNames)));
    col = @(p) find(startsWith(vn, p), 1);
    tt = T{:, col("time")};
    if ~isdatetime(tt), tt = datetime(string(tt), 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); end
    if isempty(tt.TimeZone), tt.TimeZone = 'UTC'; end
    L = table(tt, T{:, col("hm0")}, T{:, col("tm01")}, T{:, col("meandir")}, ...
        'VariableNames', {'TimeUTC','Hm0','Tm01','Meandir'});
    L = sortrows(L(~isnat(L.TimeUTC), :), 'TimeUTC');
end

function C = crosstab_table(a, b)
    ga = ["ice","miz","open"];  gb = ["ice","miz","open","unknown"];
    M = zeros(numel(ga), numel(gb));
    for i = 1:numel(ga), for j = 1:numel(gb), M(i,j) = nnz(a == ga(i) & b == gb(j)); end, end
    C = array2table(M, 'RowNames', "pipeline_" + ga, 'VariableNames', "SIC_" + gb);
end

function DSP = add_era5(DSP, Results, fn, max_dt_h, max_gap_s)
% SIC at the buoy, ice-edge normal (bearing of decreasing SIC = where waves enter the ice
% from, nautical), and 10 m wind (coming from). Validation only.
    n = height(DSP);
    DSP.SIC = nan(n,1);  DSP.ice_normal_from_deg = nan(n,1);
    DSP.wind_from_deg = nan(n,1);  DSP.wind_speed = nan(n,1);
    if ~isfile(fn), warning('%s not found: ERA5 validation skipped.', fn); return; end
    info = ncinfo(fn);  vn = {info.Variables.Name};
    lat = double(ncread(fn, pickvar(vn, {'latitude','lat'})));
    lon = double(ncread(fn, pickvar(vn, {'longitude','lon'})));
    t   = nc_time(fn, pickvar(vn, {'valid_time','time'}));
    sz  = [numel(lon), numel(lat), numel(t)];
    if ~ismember('siconc', vn), warning('siconc not in %s.', fn); return; end
    SIC = orient3(double(ncread(fn, 'siconc')), sz);
    U = []; V = [];
    if all(ismember({'u10','v10'}, vn))
        U = orient3(double(ncread(fn, 'u10')), sz);  V = orient3(double(ncread(fn, 'v10')), sz);
    end
    vnR = Results.Properties.VariableNames;
    latc = vnR(ismember(lower(vnR), {'lat','latitude','lat_deg','lat_gnss'}));
    lonc = vnR(ismember(lower(vnR), {'lon','longitude','lon_deg','lon_gnss'}));
    if isempty(latc) || isempty(lonc), warning('Buoy lat/lon not found in Results: ERA5 skipped.'); return; end
    tR = posixtime(Results.TimeUTC);  la = Results.(latc{1});  lo = Results.(lonc{1});
    fprintf('ERA5: %s to %s, %d x %d grid; buoy lat/lon from Results.%s/%s\n', t(1), t(end), ...
        sz(1), sz(2), latc{1}, lonc{1});
    miss = struct('gnss', 0, 'time', 0, 'nan', 0);
    for i = 1:n
        [g, k] = min(abs(tR - posixtime(DSP.BurstMid(i))));
        if g > max_gap_s || isnan(la(k)), miss.gnss = miss.gnss + 1; continue; end
        blat = la(k);  blon = lo(k);  if max(lon) > 180, blon = mod(blon, 360); end
        [~, ia] = min(abs(lat - blat));  [~, io] = min(abs(lon - blon));
        [dtm, it] = min(abs(t - DSP.BurstMid(i)));
        if hours(dtm) > max_dt_h, miss.time = miss.time + 1; continue; end
        DSP.SIC(i) = SIC(io, ia, it);
        if isnan(DSP.SIC(i)), miss.nan = miss.nan + 1; end     % NaN = land/missing in ERA5
        % gradient by central differences (per km)
        io2 = min(max(io + [-1 1], 1), sz(1));  ia2 = min(max(ia + [-1 1], 1), sz(2));
        s_e = SIC(io2(2), ia, it) - SIC(io2(1), ia, it);
        s_n = SIC(io, ia2(2), it) - SIC(io, ia2(1), it);
        dx = 111.32 * cosd(blat) * abs(lon(io2(2)) - lon(io2(1)));
        dy = 110.57 * abs(lat(ia2(2)) - lat(ia2(1))) * sign(lat(ia2(2)) - lat(ia2(1)));
        ge = s_e / max(dx, eps);  gn = s_n / dy;
        if all(isfinite([ge gn])) && hypot(ge, gn) > 0
            DSP.ice_normal_from_deg(i) = mod(atan2d(-ge, -gn), 360);   % toward lower SIC
        end
        if ~isempty(U)
            u = U(io, ia, it);  v = V(io, ia, it);
            DSP.wind_speed(i) = hypot(u, v);
            DSP.wind_from_deg(i) = mod(atan2d(-u, -v), 360);           % meteorological, coming from
        end
    end
    fprintf(['ERA5 SIC found for %d of %d bursts (no GNSS position within gap: %d, ERA5 time > %g h ' ...
        'away: %d, NaN cell: %d)\n'], nnz(~isnan(DSP.SIC)), n, miss.gnss, max_dt_h, miss.time, miss.nan);
end

function name = pickvar(vn, cands)
    k = find(ismember(cands, vn), 1);
    if isempty(k), error('None of {%s} found in the ERA5 file.', strjoin(cands, ', ')); end
    name = cands{k};
end

function t = nc_time(fn, name)
    v = double(ncread(fn, name));  v = v(:);
    try, u = ncreadatt(fn, name, 'units'); catch, u = 'seconds since 1970-01-01'; end
    tok = regexp(u, '^\s*(\w+)\s+since\s+(\d{4}-\d{1,2}-\d{1,2})', 'tokens', 'once');
    t0 = datetime(tok{2}, 'InputFormat', 'yyyy-M-d', 'TimeZone', 'UTC');
    switch lower(tok{1})
        case {'seconds','second','s'}, t = t0 + seconds(v);
        case {'minutes','minute'},     t = t0 + minutes(v);
        case {'hours','hour','h'},     t = t0 + hours(v);
        case {'days','day'},           t = t0 + days(v);
        otherwise, error('Unrecognised time units "%s".', u);
    end
end

function A = orient3(A, n)
    A = squeeze(A);  sz = size(A);  sz(end+1:3) = 1;
    if isequal(sz(1:3), n), return; end
    P = perms(1:3);
    for k = 1:size(P, 1)
        if isequal(sz(P(k,:)), n), A = permute(A, P(k,:)); return; end
    end
    error('ERA5 field size [%s] does not match lon/lat/time [%s].', num2str(sz), num2str(n));
end

function write_summary(S, fn, sheet)
% Append S to the summary sheet, aligning columns by name (union of old and new columns).
    if isfile(fn) && any(sheetnames(fn) == sheet)
        old = readtable(fn, 'Sheet', sheet, 'TextType', 'string');
        for v = setdiff(string(S.Properties.VariableNames), string(old.Properties.VariableNames))
            old.(v) = missing_like(S.(v), height(old));
        end
        for v = setdiff(string(old.Properties.VariableNames), string(S.Properties.VariableNames))
            S.(v) = missing_like(old.(v), 1);
        end
        S = S(:, old.Properties.VariableNames);
        for v = string(S.Properties.VariableNames)              % harmonise types for vertcat
            if isnumeric(old.(v)) && ~isnumeric(S.(v)), S.(v) = double(string(S.(v))); end
            if (isstring(old.(v)) || iscellstr(old.(v))) && ~isstring(S.(v)), S.(v) = string(S.(v)); end
            if isstring(S.(v)) && ~isstring(old.(v)), old.(v) = string(old.(v)); end
            if isdatetime(S.(v)) && ~isdatetime(old.(v)), old.(v) = datetime(string(old.(v))); end
        end
        S = [old; S];
    end
    writetable(S, fn, 'Sheet', sheet, 'WriteMode', 'overwritesheet');
end

function m = missing_like(x, n)
    if isnumeric(x) || islogical(x), m = nan(n, 1);
    elseif isdatetime(x),           m = NaT(n, 1);
    else,                           m = strings(n, 1);
    end
end
