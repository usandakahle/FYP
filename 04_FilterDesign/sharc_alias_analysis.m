%% sharc_alias_analysis.m
% Justifies the anti-aliasing filter in run_sharc_jack.m (V1.3) from SB04 data
% and produces figure:
%    stopband: attenuation required by the data vs attenuation achieved
%    passband: gain across the wave band for the 60 dB and 40 dB designs

% Pipeline: Kaiser FIR (design_aa), pass 0.5 Hz, stop fs_dec - 0.5 = 2.625 Hz,
% applied with filtfilt (power gain |H|^4), then downsample by M = 32 (100 -> 3.125 Hz).

% Required attenuation at each source frequency f_src that folds onto a wave-band
% frequency f_d:   A_req = 10log10( S(f_src) / S(f_d) ) + margin,
% both PSDs from the raw 100 Hz record (in which the wave band is not aliased).

% INPUT: raw-count packets (Heave ~ 8192 at rest), not the 2 Hz low-passed files.
% Requires the Signal Processing Toolbox.

clear; close all; clc;

%% SETTINGS
cfg.dataDir   = fullfile('SB04','sb04_packets');
cfg.filePat   = '*.csv';
cfg.fs        = 100;            % Hz
cfg.M         = 32;             % opts.M
cfg.f_pass    = 0.5;            % opts.f_max_valid
cfg.atten     = [60 40];        % [pipeline, comparison] stopband spec (dB, per pass)
cfg.accScale  = 9.81/8192;      % counts -> m/s^2
cfg.gyrScale  = 1/65.5;         % counts -> deg/s
cfg.g         = 9.81;
cfg.channels  = {'Surge','Sway','Heave','Roll','Pitch'};
cfg.isAccel   = [true true true false false];
cfg.waveBand  = [0.03 0.5];     % opts.f_min .. opts.f_max_valid
cfg.Nfft_dec  = 512;            % opts.Nfft (df = 0.0061 Hz)
cfg.margin_dB = 10;
cfg.pctPacket = 100;            % 100 = worst packet, 95 = 95th percentile
cfg.minLen    = 0.9*180000;
cfg.outDir    = 'alias_analysis_out';

if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end
fs_dec = cfg.fs/cfg.M;
L      = [cfg.Nfft_dec*cfg.M, cfg.Nfft_dec];
gainFF = @(h, f) abs(freqz(h, 1, f, cfg.fs)).^4;      % filtfilt power gain

h = cell(1,2);
for i = 1:2, h{i} = design_aa(cfg.fs, fs_dec, cfg.f_pass, cfg.atten(i)); end

%% LOAD PACKETS
files = dir(fullfile(cfg.dataDir, cfg.filePat));
assert(~isempty(files), 'No files found in %s', cfg.dataDir);
packets = struct('x',{});
for iF = 1:numel(files)
    fn   = fullfile(files(iF).folder, files(iF).name);
    opts = detectImportOptions(fn, 'VariableNamingRule','preserve');
    opts.SelectedVariableNames = [{'timeSeriesCount'}, cfg.channels];
    T    = readtable(fn, opts);
    for id = unique(T.timeSeriesCount).'
        rows = T.timeSeriesCount == id;
        if nnz(rows) < cfg.minLen, continue; end
        for c = 1:numel(cfg.channels)
            raw = double(T.(cfg.channels{c})(rows));  raw = raw(2:end);
            if cfg.isAccel(c)
                x = raw*cfg.accScale;
                if strcmp(cfg.channels{c},'Heave'), x = x - cfg.g; end
            else
                x = raw*cfg.gyrScale;
            end
            p.x.(cfg.channels{c}) = detrend(x, 2);
        end
        packets(end+1) = p; %#ok<SAGROW>
    end
end
nP = numel(packets);

%% ALIAS MAP
fgrid = (0:L(2)/2).' * fs_dec/L(2);
fd    = fgrid(fgrid >= cfg.waveBand(1) & fgrid <= cfg.waveBand(2));
kmax  = floor((cfg.fs/2 + cfg.waveBand(2))/fs_dec);
src = []; dst = [];
for k = 1:kmax
    for s = [-1 1]
        fsrc = k*fs_dec + s*fd;  ok = fsrc > 0 & fsrc <= cfg.fs/2;
        src = [src; fsrc(ok)];  dst = [dst; find(ok)]; %#ok<AGROW>
    end
end
[src, o] = sort(src);  dst = dst(o);
Aach_src = -10*log10(gainFF(h{1}, src));

%% PROCESS
w = @(n) hann(n,'periodic');
Areq = zeros(numel(src), nP, numel(cfg.channels));
aliasRatio = zeros(nP, numel(cfg.channels));
measExcess = zeros(nP, numel(cfg.channels));
Hpass = gainFF(h{1}, fd);
for c = 1:numel(cfg.channels)
    for iP = 1:nP
        x0 = packets(iP).x.(cfg.channels{c});
        [P0, f0] = pwelch(x0, w(L(1)), L(1)/2, L(1), cfg.fs);
        x2 = downsample(filtfilt(h{1}, 1, x0), cfg.M);
        [P2, f2] = pwelch(x2, w(L(2)), L(2)/2, L(2), fs_dec);
        Strue = interp1(f0, P0, fd);
        Ssrc  = interp1(f0, P0, src);
        Areq(:,iP,c)    = 10*log10(Ssrc ./ Strue(dst)) + cfg.margin_dB;
        Salias          = accumarray(dst, gainFF(h{1}, src).*Ssrc, [numel(fd) 1]);
        aliasRatio(iP,c) = max(10*log10(Salias ./ Strue));
        measExcess(iP,c) = max(abs(10*log10(interp1(f2, P2, fd) ./ (Strue .* Hpass))));
    end
end

% Envelope of the requirement: chosen packet percentile, worst channel.
Areq_env = -Inf(numel(src),1);  chWorst = zeros(numel(src),1);
for c = 1:numel(cfg.channels)
    a = pct(Areq(:,:,c), cfg.pctPacket);
    up = a > Areq_env;  Areq_env(up) = a(up);  chWorst(up) = c;
end
[AreqMax, iMax] = max(Areq_env);

%% FIGURE
ff  = linspace(0, cfg.fs/2, 20001).';
A   = [-10*log10(gainFF(h{1}, ff)), -10*log10(gainFF(h{2}, ff))];
fpb = linspace(0, 0.6, 601).';
G   = [10*log10(gainFF(h{1}, fpb)), 10*log10(gainFF(h{2}, fpb))];

% Break the requirement line between alias bands.
gap = [false; diff(src) > 2*(fs_dec/L(2))];
xs = src;  ys = Areq_env;
xs = repelem(xs, 1 + gap);  ys = repelem(ys, 1 + gap);
ins = find(gap) + (0:nnz(gap)-1)';      % positions of inserted copies
xs(ins) = NaN;  ys(ins) = NaN;

if cfg.pctPacket == 100, pk = 'worst packet'; else, pk = sprintf('%g^{th} pct packet', cfg.pctPacket); end
fig = figure('Units','centimeters','Position',[2 2 16 13],'Color','w');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing','compact','Padding','compact');

% (a) Stopband
ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
for k = 1:kmax
    a = k*fs_dec - cfg.f_pass;  b = min(k*fs_dec + cfg.f_pass, cfg.fs/2);
    patch(ax1, [a b b a], [0 0 200 200], [0.92 0.92 0.92], 'EdgeColor','none', 'HandleVisibility','off');
end
plot(ax1, ff, A(:,1), 'k-',  'LineWidth', 1.0, 'DisplayName', sprintf('Achieved, %g dB design', cfg.atten(1)));
plot(ax1, ff, A(:,2), 'k--', 'LineWidth', 0.8, 'DisplayName', sprintf('Achieved, %g dB design', cfg.atten(2)));
plot(ax1, xs, ys, '.-', 'Color',[0 0.35 0.75], 'MarkerSize', 6, ...
     'DisplayName', sprintf('Required (%s, %g dB margin)', pk, cfg.margin_dB));
text(ax1, src(iMax), AreqMax + 8, sprintf('max %.0f dB', AreqMax), ...
     'HorizontalAlignment','center', 'FontSize', 8, 'Color', [0 0.35 0.75]);
xlim(ax1, [0 cfg.fs/2]);  ylim(ax1, [0 160]);
xlabel(ax1, 'Frequency before decimation (Hz)');  ylabel(ax1, 'Attenuation (dB)');
title(ax1, '(a) Bands folding into 0.03-0.5 Hz (shaded): required vs achieved', 'FontWeight','normal');
legend(ax1, 'Location','southeast', 'FontSize', 8);

% (b) Passband
ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
patch(ax2, [cfg.waveBand fliplr(cfg.waveBand)], [-1 -1 1 1], [0.92 0.92 0.92], ...
      'EdgeColor','none', 'HandleVisibility','off');
yline(ax2, 0, 'k:', 'HandleVisibility','off');
plot(ax2, fpb, G(:,1), 'k-',  'LineWidth', 1.0, 'DisplayName', sprintf('%g dB design', cfg.atten(1)));
plot(ax2, fpb, G(:,2), 'k--', 'LineWidth', 0.8, 'DisplayName', sprintf('%g dB design', cfg.atten(2)));
xlim(ax2, [0 0.6]);  ylim(ax2, [-0.4 0.4]);
xlabel(ax2, 'Frequency (Hz)');  ylabel(ax2, 'Power gain (dB)');
title(ax2, '(b) Wave band (shaded): passband ripple', 'FontWeight','normal');
legend(ax2, 'Location','southwest', 'FontSize', 8);

set([ax1 ax2], 'FontSize', 9);
exportgraphics(fig, fullfile(cfg.outDir, 'alias_justification.pdf'), 'ContentType','vector');
exportgraphics(fig, fullfile(cfg.outDir, 'alias_justification.png'), 'Resolution', 300);

%% CONSOLE DEBUG
fprintf('\n%d packets, %d channels, %s, %g dB margin\n', nP, numel(cfg.channels), strrep(pk,'^{th}','th'), cfg.margin_dB);
fprintf('Max required attenuation: %.1f dB at %.2f Hz (%s)\n', AreqMax, src(iMax), cfg.channels{chWorst(iMax)});
for i = 1:2
    Aa  = -10*log10(gainFF(h{i}, src));
    sur = min(Aa - Areq_env);
    pbl = max(abs(10*log10(gainFF(h{i}, fd))));
    fprintf('%g dB design: %3d taps | min surplus %.1f dB | max passband error %.3f dB (%.2f%% energy)\n', ...
        cfg.atten(i), numel(h{i}), sur, pbl, 100*(10^(pbl/10)-1));
end
fprintf('Current design, all channels: aliased/true <= %.1f dB, measured excess <= %.1e dB\n', ...
    max(aliasRatio(:)), max(measExcess(:)));

%% LOCAL FUNCTIONS
function h = design_aa(fs, fs_dec, fpass, atten_db)
% Copied from run_sharc_jack.m.
    fstop = fs_dec - fpass;
    [nk, Wn, beta, ftype] = kaiserord([fpass fstop], [1 0], [0.01 10^(-atten_db/20)], fs);
    nk = nk + mod(nk, 2);
    h = fir1(nk, Wn, ftype, kaiser(nk+1, beta), 'noscale');
end

function q = pct(X, p)
% Percentile along dim 2.
    X = sort(X, 2);  n = size(X,2);
    q = X(:, max(1, min(n, round(p/100*n))));
end