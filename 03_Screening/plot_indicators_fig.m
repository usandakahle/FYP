%% plot_indicators.m  
% platform indicators for SB04 and SB06
% Reads <BUOY>_diagnostics.csv (written by sharc_explore.m) and plots, for each buoy:
% maximum yaw rate        reference 5 deg/s  (rotation onset rule)
% kurtosis, a_z > 2 Hz    reference 3 (Gaussian) and 3.5 (impact limit)
% tilt gain |accel/gyro|  reference 0.5 (rigid floe above)
% spectral tail exponent  reference band n = 4-5 (open water sea)
% SB04 blue, SB06 green. PNG only, written to OUTDIR as <buoy>_indicators.png.

clc;
%% SETTINGS
BUOYS   = ["SB04", "SB06"];
CSV     = @(b) b + "_diagnostics.csv";     
COLS    = struct('SB04', [0.00 0.45 0.74], 'SB06', [0.20 0.60 0.30]);
OUTDIR  = "IndicatorFigures";
DPI     = 200;
LW      = 1.0;                                % line width
FS      = 10;                                 % font size
TICK_H  = 6;                                  % hours between time ticks
REF     = [0.45 0.45 0.45];                   % reference lines

YAW_REF   = 5;        % deg/s
KURT_G    = 3;        % Gaussian
KURT_LIM  = 3.5;      % impact limit
GAIN_REF  = 0.5;      % rigid floe above
TAIL_OW   = [4 5];    % open-water tail exponents

%% PLOT
set(groot, 'defaultFigureColor', 'w');
if ~exist(OUTDIR, 'dir'), mkdir(OUTDIR); end

for b = BUOYS
    fn = CSV(b);
    if ~isfile(fn), warning('%s not found: %s skipped.', fn, b); continue; end
    X = readtable(fn);
    t = X.tmid;
    if ~isdatetime(t), t = datetime(string(t), 'InputFormat', 'dd-MMM-yyyy HH:mm:ss'); end
    [t, o] = sort(t);  X = X(o, :);
    gain = mean(abs([X.tilt_gain_x X.tilt_gain_y]), 2, 'omitnan');
    c = COLS.(char(b));

    % time ticks every TICK_H hours, on whole hours
    t0 = dateshift(min(t), 'start', 'hour');  t0 = t0 - hours(mod(hour(t0), TICK_H));
    t1 = dateshift(max(t), 'end', 'hour');
    tk = t0:hours(TICK_H):t1 + hours(TICK_H);

    fig = figure('Color', 'w', 'Position', [60 60 1100 820], 'Visible', 'off');
    tl = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax = gobjects(4, 1);

    % max yaw rate
    ax(1) = nexttile(tl);  hold(ax(1), 'on');
    plot(ax(1), t, X.yaw_max, '-', 'Color', c, 'LineWidth', LW);
    yline(ax(1), YAW_REF, '--', sprintf('Threshold %g deg/s', YAW_REF), 'Color', REF, ...
        'LabelHorizontalAlignment', 'left', 'FontSize', FS - 1);
    set(ax(1), 'YScale', 'log');
    ylabel(ax(1), 'Max |Yaw Rate| (deg/s)');

    % kurtosis
    ax(2) = nexttile(tl);  hold(ax(2), 'on');
    plot(ax(2), t, X.kurt_hf, '-', 'Color', c, 'LineWidth', LW);
    yline(ax(2), KURT_G, '-', 'Gaussian (3)', 'Color', REF, 'LabelHorizontalAlignment', 'left', 'FontSize', FS - 1);
    yline(ax(2), KURT_LIM, '--', sprintf('Impact Limit %g', KURT_LIM), 'Color', REF, 'LabelHorizontalAlignment', 'right', 'FontSize', FS - 1);
    set(ax(2), 'YScale', 'log');
    ylabel(ax(2), 'Kurtosis of a_z > 2 Hz (-)');

    % tilt gain
    ax(3) = nexttile(tl);  hold(ax(3), 'on');
    plot(ax(3), t, gain, '-', 'Color', c, 'LineWidth', LW);
    ylim(ax(3), [0 1.05]);
    ylabel(ax(3), 'Tilt Gain |Accel/Gyro| (-)');

    % tail exponent
    ax(4) = nexttile(tl);  hold(ax(4), 'on');
    plot(ax(4), t, X.tail_n, '-', 'Color', c, 'LineWidth', LW);
    ylim(ax(4), [0 max(30, ceil(max(X.tail_n, [], 'omitnan')/5)*5)]);
    ylabel(ax(4), 'Tail Exponent n (C_{zz} \propto f^{-n})');

    % common time axis, grid lines, panel labels
    lab = ["(a)", "(b)", "(c)", "(d)"];
    for k = 1:4
        grid(ax(k), 'on');  box(ax(k), 'on');
        set(ax(k), 'FontSize', FS, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'GridAlpha', 0.25);
        xlim(ax(k), [min(t) - minutes(30), max(t) + minutes(30)]);
        xticks(ax(k), tk);
        if k < 4
            ax(k).XTickLabel = [];
        else
            xtickformat(ax(k), 'dd-MMM HH:mm');
            xtickangle(ax(k), 30);
            xlabel(ax(k), 'Time (UTC)');
        end
        text(ax(k), 0.005, 0.90, lab(k), 'Units', 'normalized', 'FontSize', FS);
    end
    linkaxes(ax, 'x');

    out = fullfile(OUTDIR, lower(b) + "_indicators.png");
    exportgraphics(fig, out, 'Resolution', DPI, 'BackgroundColor', 'white');
    close(fig);
    fprintf('%s: %d packets, %s to %s -> %s\n', b, height(X), ...
        string(min(t), 'dd-MMM HH:mm'), string(max(t), 'dd-MMM HH:mm'), out);
end
