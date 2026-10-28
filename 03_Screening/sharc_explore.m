%% sharc_direction_diagnostics.m
% Independent per packet diagnosis of SHARC directional data. 

%   T0 Is the data intact?          fs, gaps, saturation, flat-lines
%   T1 How is the buoy oriented?    gravity axis and sign (z up/down), static tilt
%   T2 Is the buoy rotating?        gyro bias, net yaw over the packet, max yaw rate
%   T3 Do the axes agree?           accel tilt vs gyro tilt: 
%                                   acceleration cancel, so Surge/Sway carry no slope)
%   T4 Is there directional info?   heave-slope coherence r1^2 vs its 95 % noise level,
%                                   co/quad ratio (slopes must be 90 deg from heave),
%                                   check ratio R = (Cxx+Cyy)/(k^2 Czz) (~1)
%   T5 Which direction, in which frame?
%                                   body-frame direction at the PEAK band from gyro
%                                   slopes and from accel slopes; Czz-weighted band
%                                   mean (as the pipeline) and the share of weight
%                                   below 0.05 Hz; earth frame with GNSS yaw, for both
%                                   handedness conventions
%   T6 Is the spectrum trustworthy? Hs, fp, SNR at the peak vs the measured accel noise
%                                   floor (593 ug/rtHz), two-system split
%   T7 Regime features              tail exponent, hf fraction, kurtosis > 2 Hz,
%                                   tilt transfer ratio (T3)
% Output: one row per packet (CSV), a verdict per packet, 7 summary figures and a
% printed decision summary. Only needs the packet CSVs (+ optional GNSS yaw / reference).

clear; clc; close all;

%% SETTINGS
BUOY_ID    = "SB04";                     % change for Sb06
PACKET_DIR = BUOY_ID + "/" + lower(BUOY_ID) + "_packets";   % packet_*.csv
SCALE_ACC  = 9.80665/8192;               % counts -> m/s^2 (+/-4 g)
SCALE_GYR  = 1/65.5;                     % counts -> deg/s (+/-500 dps)
YAW_FILE   = BUOY_ID + "_yaw_state.mat"; % from build_yaw_checkpoint 
YAW_MAX_GAP_MIN = 20;
WINDOW     = [datetime(2022,7,1,'TimeZone','UTC'), datetime(2022,8,1,'TimeZone','UTC')];
NA         = 593e-6 * 9.80665;           % accel noise density, m/s^2/rtHz (6 h test)
BAND       = [0.03 0.5];                 % wave band, Hz
TILT_BAND  = [0.03 0.15];                % band for the accel-vs-gyro tilt test
NFFT       = 512;                        % at ~3.2 Hz: ~158 s segments, 50 % overlap
OUT        = BUOY_ID + "_diagnostics";
g          = 9.80665;

%% FILES 
files = dir(fullfile(PACKET_DIR, "packet_*.csv"));
if isempty(files), error('No packet_*.csv in %s', PACKET_DIR); end
tstart = NaT(numel(files),1,'TimeZone','UTC');
for k = 1:numel(files)
    tok = regexp(files(k).name, 'packet_(\d+)_(.+)\.csv', 'tokens', 'once');
    if ~isempty(tok), tstart(k) = datetime(tok{2}, 'InputFormat','yyyy-MM-dd_HH-mm-ss', 'TimeZone','UTC'); end
end
keep = tstart >= WINDOW(1) & tstart <= WINDOW(2);
files = files(keep);  [~, o] = sort(tstart(keep));  files = files(o);
fprintf('%d packets in window.\n', numel(files));

% optional GNSS yaw
G = [];
if strlength(YAW_FILE) > 0 && isfile(YAW_FILE)
    S = load(YAW_FILE, 'Results');
    G = table(S.Results.TimeUTC, S.Results.YawEstimate_deg, 'VariableNames', {'t','yaw'});
    if isempty(G.t.TimeZone), G.t.TimeZone = 'UTC'; end
    G = G(~isnan(G.yaw) & ~isnat(G.t), :);
end

%% PER PACKET 
rows = struct([]);
for k = 1:numel(files)
    fp_ = fullfile(files(k).folder, files(k).name);
    try
        T = readtable(fp_);
        if ~isdatetime(T.DateTimeUTC)
            T.DateTimeUTC = datetime(T.DateTimeUTC, 'InputFormat','dd-MMM-yyyy HH:mm:ss', 'TimeZone','UTC');
        elseif isempty(T.DateTimeUTC.TimeZone)
            T.DateTimeUTC.TimeZone = 'UTC';
        end
        D = diag_packet(T, SCALE_ACC, SCALE_GYR, NA, BAND, TILT_BAND, NFFT, g);
    catch ME
        warning('%s: %s', files(k).name, ME.message);  continue;
    end
    D.file = string(files(k).name);
    % heading from GNSS yaw (nearest fix)
    D.yaw = NaN;
    if ~isempty(G)
        [gap, iy] = min(abs(G.t - D.tmid));
        if minutes(gap) <= YAW_MAX_GAP_MIN, D.yaw = mod(G.yaw(iy), 360); end
    end
    D.verdict = verdict(D);
    % fields needed only for the verdict, not written to the CSV
    D = rmfield(D, {'n_gaps','sat_acc','sat_gyr','flat_s','snr_peak','r1sq_thr','n_peaks','dir_low','dir_high'});
    if isempty(rows), rows = D; else, rows(end+1) = D; end %#ok<AGROW>
    fprintf('[%3d/%d] %s  %-22s r1^2 %.2f  tilt gain %+.2f\n', k, numel(files), ...
        datestr(D.tmid, 'dd HH:MM'), D.verdict, D.r1sq, D.tilt_gain_x);
end
if isempty(rows), error('No packets processed.'); end
X = struct2table(rows);
X = movevars(X, {'file','tmid','verdict'}, 'Before', 1);
Xo = X;  Xo.tmid.TimeZone = '';
writetable(Xo, OUT + ".csv");
fprintf('\nWrote %s.csv (%d packets)\n', OUT, height(X));

%%  FIGURES 
t = X.tmid;
f1 = figure('Name','T0-T1 data and orientation');
plotm(t, X.fs_ts, 'fs from timestamps (Hz)', [100]);

f2 = figure('Name','T2 rotation');
tiledlayout(f2, 2, 1);
nexttile; plotm(t, X.yaw_net, 'net yaw over packet (deg)', [-90 0 90]);
nexttile; plotm(t, X.yaw_max, 'max |yaw rate| (deg/s)', 30);

f3 = figure('Name','T3 accel vs gyro tilt');
plotm(t, [X.tilt_gain_x X.tilt_gain_y], 'tilt transfer gain (accel/gyro)', [0 1]);
legend('x', 'y'); title('1 = rigid platform (ice), 0 = wave-following (open water), < 0 = sign/mapping error');

f4 = figure('Name','T4 directional information');
tiledlayout(f4, 3, 1);
nexttile; plotm(t, X.r1sq, 'heave-slope coherence r_1^2', []);
nexttile; plotm(t, X.coquad, 'co/quad ratio (want << 1)', [0.5]);
nexttile; plotm(t, X.check_R, 'check ratio R (want ~1)', [0.5 1 2]); set(gca, 'YScale', 'log');

f7 = figure('Name','T6-T7 spectrum and regime features');
tiledlayout(f7, 2, 2);
nexttile; plotm(t, [X.Hs X.Hs_noise], 'Hs (m)', []); legend('Hs', 'noise-floor Hs');
nexttile; plotm(t, X.tail_n, 'tail exponent n', [4 5 7]);
nexttile; plotm(t, X.kurt_hf, 'kurtosis > 2 Hz', 3); set(gca,'YScale','log');
nexttile; plotm(t, X.w_below_005, 'Czz weight below 0.05 Hz', 0.3);

for h = [f1 f2 f3 f4 f7]
    exportgraphics(h, OUT + "_" + regexprep(h.Name, '[^A-Za-z0-9]+', '_') + ".png", 'Resolution', 150);
end

%%  DECISION SUMMARY 
fprintf('\n DECISION SUMMARY (%s, %d packets) \n', BUOY_ID, height(X));
fprintf('Verdicts:\n');  disp(groupcounts(X, 'verdict'));
fprintf('T4 co/quad > 0.5 (phase problem) in %d packets\n', nnz(X.coquad > 0.5));
fprintf('=====================================================================\n');

%% LOCAL FUNCTIONS
function D = diag_packet(T, sa, sg, NA, BAND, TILT_BAND, NFFT, g)
    t = T.DateTimeUTC;  n = height(T);
    D = struct();
    D.tmid = t(1) + (t(end) - t(1))/2;
    span = seconds(t(end) - t(1));
    D.n = n;
    D.fs_ts = (n - 1) / max(span, eps);
    fs = D.fs_ts;  if abs(fs/100 - 1) > 0.1, fs = 100; end
    dtt = seconds(diff(t));  D.n_gaps = nnz(dtt > 1.5);
    raw = double([T.Surge T.Sway T.Heave T.Roll T.Pitch T.Yaw]);
    raw = fillmissing(raw, 'linear');
    D.sat_acc = mean(any(abs(raw(:,1:3)) >= 32700, 2));
    D.sat_gyr = mean(any(abs(raw(:,4:6)) >= 32700, 2));
    D.flat_s  = max_run(raw(:,3)) / fs;
    acc = sa * raw(:,1:3);  gyr = sg * raw(:,4:6);           % m/s^2, deg/s
   
    % T1 orientation
    gm = mean(acc, 1);
    zs = sign(gm(3));
    
    % T2 rotation
    D.yaw_net = sum(gyr(:,3) - median(gyr(:,3))) / fs;     % net rotation beyond a constant bias
    D.yaw_max = max(abs(movmedian(gyr(:,3) - median(gyr(:,3)), 5)));
    
    % kurtosis of heave accel above 2 Hz (collisions)
    [b, a] = butter(4, 2/(fs/2), 'high');
    hf = filtfilt(b, a, acc(:,3) - mean(acc(:,3)));
    D.kurt_hf = mean(hf.^4) / max(mean(hf.^2)^2, eps);
    
    % decimate to ~fs/32 (two stages)
    dec = @(x) decimate(decimate(x, 4, 'fir'), 8, 'fir');
    Y = [dec(acc(:,3)) dec(acc(:,1)) dec(acc(:,2)) dec(deg2rad(gyr(:,1))) dec(deg2rad(gyr(:,2))) dec(deg2rad(gyr(:,3)))];
    fsd = fs / 32;
    [f, Gm, nseg] = welch_matrix(Y, fsd, NFFT);
    w = 2*pi*f;  w(1) = Inf;                                % DC bin - all factors 0 (no NaN)
   
    % channel factors: heave disp, accel tilt, gyro angle (physical signs from z orientation)
    if zs == 0, zs = 1; end
    hZ = -1 ./ w.^2;  hA = ones(size(f)) / g;  hG = -1i ./ w;   % 1/(i w), finite at w = Inf
    sx_pitch = -zs;  sy_roll = zs;                         % z up: Sx = -pitch, Sy = +roll
    
    % spectral entries (1 az, 2 ax, 3 ay, 4 roll(x), 5 pitch(y), 6 yaw)
    Gt = @(i, j, hi, hj) conj(hi) .* hj .* Gm(:, i, j);
    Czz = real(Gt(1,1,hZ,hZ));
    Sgx = {5, sx_pitch*hG};  Sgy = {4, sy_roll*hG};        % gyro slopes
    Sax = {2, hA};           Say = {3, hA};                % accel slopes
    GZ = @(c) Gt(1, c{1}, hZ, c{2});
    GG = @(c, d) Gt(c{1}, d{1}, c{2}, d{2});
    inb = f >= BAND(1) & f <= BAND(2);
    df = f(2) - f(1);
    
    % T6 spectrum
    D.Hs = 4*sqrt(sum(Czz(inb)) * df);
    En = NA^2 ./ w.^4;
    D.Hs_noise = 4*sqrt(sum(En(inb)) * df);
    sel = f >= 0.04 & f <= 0.3;
    [~, kp] = max(Czz .* sel);  fp = f(kp);  D.fp = fp;
    D.snr_peak = Czz(kp) / En(kp);
    pk = abs(f - fp) <= 0.25*fp & inb;                     % peak band
    
    % T3 accel tilt vs gyro tilt
    tb = f >= TILT_BAND(1) & f <= TILT_BAND(2);
    D.tilt_gain_x = real(sum(GG(Sgx, Sax) .* tb) / sum(real(GG(Sgx, Sgx)) .* tb));
    D.tilt_gain_y = real(sum(GG(Sgy, Say) .* tb) / sum(real(GG(Sgy, Sgy)) .* tb));
    
    % T4 directional information (peak band)
    Qx = sum(GZ(Sgx) .* pk);  Qy = sum(GZ(Sgy) .* pk);
    Cxx = sum(real(GG(Sgx, Sgx)) .* pk);  Cyy = sum(real(GG(Sgy, Sgy)) .* pk);
    Cz = sum(Czz .* pk);
    a1 = -imag(Qx) / sqrt(Cz*(Cxx + Cyy));  b1 = -imag(Qy) / sqrt(Cz*(Cxx + Cyy));
    D.r1sq = a1^2 + b1^2;
    D.r1sq_thr = 1 - 0.05^(1/max(nseg - 1, 1));
    D.coquad = hypot(real(Qx), real(Qy)) / max(hypot(imag(Qx), imag(Qy)), eps);
    k = (w.^2) / g;  k(1) = 0;
    D.check_R = (Cxx + Cyy) / sum(k.^2 .* Czz .* pk);
    
    % low-frequency share of the Czz weight
    wgt = Czz .* (inb & f <= 0.3);
    D.w_below_005 = sum(wgt .* (f < 0.05), 'omitnan') / sum(wgt, 'omitnan');
    
    % two systems: split at the valley between the two largest peaks
    c = movmean(Czz, 3);  c(~sel) = 0;
    [pkv, loc] = findpeaks(c, 'MinPeakProminence', 0.1*max(c));
    D.n_peaks = numel(pkv);  D.dir_low = NaN;  D.dir_high = NaN;
    if numel(pkv) >= 2
        [~, o] = sort(pkv, 'descend');  two = sort(loc(o(1:2)));
        [~, kv] = min(c(two(1):two(2)));  kv = two(1) + kv - 1;
        bands = {sel & f < f(kv), sel & f >= f(kv)};
        dd = nan(1,2);
        for q = 1:2
            m = bands{q};
            dd(q) = rel_from(atan2d(-imag(sum(GZ(Sgy) .* m)), -imag(sum(GZ(Sgx) .* m))), zs);
        end
        D.dir_low = dd(1);  D.dir_high = dd(2);
    end
    
    % T7 regime features
    floor_ = median(Czz(f >= 0.35 & f <= 0.5));
    hi = f > 1.3*fp & f <= 0.3 & Czz > 10*max(floor_, eps);
    D.tail_n = NaN;
    if nnz(hi) >= 4, pp = polyfit(log(f(hi)), log(Czz(hi)), 1);  D.tail_n = -pp(1); end
end

function [rel, rel_alt] = rel_from(phi, zs)
% phi: direction waves travel TOWARD, angle from body x toward body y (deg).
% Returns the direction waves come FROM as a clockwise bearing from the bow.
% z up (right-handed): y = port, phi is anticlockwise -> bearing = -phi.
% z down: y = starboard, phi is clockwise -> bearing = +phi.
    if zs > 0, rel = mod(180 - phi, 360);  rel_alt = mod(180 + phi, 360);
    else,      rel = mod(180 + phi, 360);  rel_alt = mod(180 - phi, 360);
    end
end

function [f, G, nseg] = welch_matrix(X, fs, nfft)
% Welch cross-spectral matrix, Hann, 50 % overlap. G(:,i,j) = mean(conj(Xi).*Xj), one-sided.
    N = size(X, 1);  C = size(X, 2);
    w = hann(nfft, 'periodic');  hop = nfft/2;
    st = 1:hop:(N - nfft + 1);  nseg = numel(st);
    if nseg < 2, error('packet too short for nfft = %d', nfft); end
    nf = nfft/2 + 1;  G = zeros(nf, C, C);
    for s = st
        seg = X(s:s+nfft-1, :);  seg = seg - mean(seg, 1);
        F = fft(seg .* w, nfft);  F = F(1:nf, :);
        for i = 1:C, for j = 1:C, G(:, i, j) = G(:, i, j) + conj(F(:, i)) .* F(:, j); end, end
    end
    G = G * 2 / (fs * sum(w.^2) * nseg);
    f = (0:nf-1)' * fs / nfft;
end

function v = verdict(D)
% First applicable explanation for a packet (ordered from data to physics).
    if D.n_gaps > 0 || D.sat_acc > 1e-4 || D.sat_gyr > 1e-4 || D.flat_s > 5 || abs(D.fs_ts/100 - 1) > 0.1
        v = "DATA: gaps/saturation";
    elseif abs(D.yaw_net) > 90 || D.yaw_max > 30
        v = "SPIN: buoy rotated";
    elseif D.kurt_hf > 100
        v = "IMPACTS";
    elseif D.snr_peak < 10
        v = "LOW SNR";
    elseif D.r1sq < D.r1sq_thr
        v = "NO DIR INFO";
    elseif D.coquad > 0.5
        v = "PHASE: co/quad";
    elseif D.check_R < 0.5 || D.check_R > 2
        v = "TILT != SLOPE";
    elseif D.w_below_005 > 0.3
        v = "LOW-F WEIGHT";
    elseif D.n_peaks >= 2 && abs(wrap180(D.dir_low - D.dir_high)) > 45
        v = "TWO SYSTEMS";
    else
        v = "OK";
    end
end

function plotm(t, y, ylab, refs)
    plot(t, y, '-o', 'MarkerSize', 3);  hold on;  grid on;  ylabel(ylab);
    for r = refs(:)', yline(r, '--', 'Color', [0.4 0.4 0.4]); end
end

function n = max_run(x)
    d = [true; diff(x) ~= 0; true];  n = max(diff(find(d)));
end

function y = wrap180(x), y = mod(x + 180, 360) - 180; end
