%% build_yaw_checkpoint.m
% Builds <BUOY_ID>_yaw_state.mat for batchBias.m.

% "YAW": course over ground (COG) between consecutive valid GNSS fixes,
% smoothed by a causal, distance weighted circular filter. It is the DRIFT BEARING,
% used as a heading proxy under the assumption that the hull weathervanes into its
% drift. That assumption fails when the buoy is frozen into a floe (no weathervaning)
% and when ice drift is set by wind/current rather than the hull.

% Also compares the buoy's COG with a reference drifter (REF_FILE, e.g. LP8) to check
% the drift-alignment assumption. Validation only; nothing from REF enters the yaw.

clear; clc;

%% SETTINGS
BUOY_ID   = "SB04";
GNSS_FILE = BUOY_ID + "/" + BUOY_ID + "GNSSTable.csv";   % DateTimeUTC, Latitude, Longitude, FixType, HDOP_*
REF_ID    = "LP8";                                        % reference drifter ("" to skip)
REF_FILE  = "RawData/" + REF_ID + ".csv";                 % 'Time (UTC)', Latitude, Longitude
WINDOW    = [NaT NaT];                                    % [start end] UTC, or NaT = whole file
MIN_FIX   = 3;            % FixType >= 3
MAX_HDOP  = 5;
MIN_DIST_M = 10;          % shorter intervals are GNSS noise
REF_MAX_GAP_S = 30*60;    % buoy-reference time match
ALPHA     = 0.3;          % yaw filter learning rate
OUT_FILE  = BUOY_ID + "_yaw_state.mat";

%%  BUOY GNSS TO COG 
op = detectImportOptions(GNSS_FILE);  op.VariableNamingRule = 'preserve';
Gt = readtable(GNSS_FILE, op);
t = Gt.DateTimeUTC;
if ~isdatetime(t), t = datetime(string(t), 'InputFormat', 'dd-MMM-yyyy HH:mm:ss'); end
t.TimeZone = 'UTC';
if ~isnat(WINDOW(1))
    w = t >= WINDOW(1) & t <= WINDOW(2);  Gt = Gt(w, :);  t = t(w);
end
[t, o] = sort(t);  Gt = Gt(o, :);
vn = Gt.Properties.VariableNames;
hdop = nan(height(Gt), 1);
if ismember('HDOP_Digit', vn), hdop = double(Gt.HDOP_Digit); end
if all(ismember({'HDOP_Digit','HDOP_Precision'}, vn))
    hdop = double(Gt.HDOP_Digit) + double(Gt.HDOP_Precision)/100;
end
lat = Gt.Latitude;  lon = Gt.Longitude;
valid = isfinite(lat) & isfinite(lon) & lat ~= 0 & lon ~= 0;
if ismember('FixType', vn), valid = valid & Gt.FixType >= MIN_FIX; end
valid = valid & ~(hdop > MAX_HDOP);
fprintf('%s: %d GNSS rows, %d valid fixes\n', BUOY_ID, height(Gt), nnz(valid));

% COG between consecutive VALID fixes (assigned to the later fix)
iv = find(valid);
[d_m, cog] = haversine_cog(lat(iv(1:end-1)), lon(iv(1:end-1)), lat(iv(2:end)), lon(iv(2:end)));
dist = nan(height(Gt), 1);  COG = nan(height(Gt), 1);  okCOG = false(height(Gt), 1);
dist(iv(2:end)) = d_m;  COG(iv(2:end)) = cog;  okCOG(iv(2:end)) = d_m >= MIN_DIST_M;
fprintf('COG intervals: %d valid of %d\n', nnz(okCOG), numel(d_m));

pkt = nan(height(Gt), 1);
if ismember('timeSeriesCount', vn), pkt = Gt.timeSeriesCount; end
Results = table(pkt, t, lat, lon, hdop, dist, COG, okCOG, 'VariableNames', ...
    {'Packet','TimeUTC','Latitude','Longitude','HDOP','Distance_m','COG_deg','ValidCOG'});

%%  YAW PROXY (causal, distance-weighted) 
refDist = median(Results.Distance_m(Results.ValidCOG), 'omitnan');
[Results.YawEstimate_deg, Results.YawConfidence] = ...
    adaptive_yaw_filter(Results.COG_deg, Results.Distance_m, Results.ValidCOG, ALPHA, refDist);
fprintf('Yaw proxy defined for %d of %d rows (ref distance %.0f m)\n', ...
    nnz(~isnan(Results.YawEstimate_deg)), height(Results), refDist);

%%  REFERENCE DRIFTER CHECK (validation only) 
Comparison_REF = table();
if strlength(REF_ID) > 0 && isfile(REF_FILE)
    op = detectImportOptions(REF_FILE);  op.VariableNamingRule = 'preserve';
    L = readtable(REF_FILE, op);
    tc = intersect({'Time (UTC)','Time_UTC_','DateTimeUTC','TimeUTC'}, L.Properties.VariableNames, 'stable');
    if isempty(tc), error('No time column in %s.', REF_FILE); end
    tl = L.(tc{1});
    if ~isdatetime(tl), tl = datetime(string(tl), 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); end
    tl.TimeZone = 'UTC';
    [tl, o] = sort(tl);  L = L(o, :);
    [dl, cl] = haversine_cog(L.Latitude(1:end-1), L.Longitude(1:end-1), L.Latitude(2:end), L.Longitude(2:end));
    tl = tl(2:end);  ok = dl >= MIN_DIST_M;  tl = tl(ok);  cl = cl(ok);
    use = find(Results.ValidCOG);
    m = false(size(use));  refc = nan(size(use));
    for j = 1:numel(use)
        [gap, k] = min(abs(seconds(tl - Results.TimeUTC(use(j)))));
        if ~isempty(gap) && gap <= REF_MAX_GAP_S, m(j) = true;  refc(j) = cl(k); end
    end
    use = use(m);  refc = refc(m);
    Comparison_REF = table(Results.TimeUTC(use), Results.COG_deg(use), refc, ...
        wrap180(Results.COG_deg(use) - refc), Results.YawEstimate_deg(use), ...
        'VariableNames', {'TimeUTC', char(BUOY_ID + "_COG_deg"), char(REF_ID + "_COG_deg"), 'COG_diff_deg', 'YawEstimate_deg'});
    z = mean(exp(1i*deg2rad(Comparison_REF.COG_diff_deg)));
    fprintf('%s vs %s drift: %d matches, mean offset %+.1f deg, circular std %.1f deg\n', BUOY_ID, REF_ID, ...
        height(Comparison_REF), rad2deg(angle(z)), rad2deg(sqrt(-2*log(max(abs(z), eps)))));
    fprintf('(similar drift says both floes move together; it does NOT show that drift = hull heading)\n');
elseif strlength(REF_ID) > 0
    warning('%s not found: reference drift check skipped.', REF_FILE);
end

%% SAVE 
meta = struct('buoy_id', BUOY_ID, 'created_utc', datetime('now', 'TimeZone', 'UTC'), ...
    'method', "GNSS course over ground, adaptive circular filter (drift-bearing heading proxy)", ...
    'alpha', ALPHA, 'ref_distance_m', refDist, 'gnss_file', GNSS_FILE, 'ref_id', REF_ID);
save(OUT_FILE, 'Results', 'Comparison_REF', 'meta');
writetable(Results, BUOY_ID + "_GNSS_COG_results.csv");
fprintf('Saved %s (buoy %s). Batch: set BUOY_ID = "%s".\n', OUT_FILE, BUOY_ID, BUOY_ID);

%% LOCAL FUNCTIONS
function [d_m, cog] = haversine_cog(lat1, lon1, lat2, lon2)
% Great-circle distance (m) and initial bearing (deg, clockwise from north).
    R = 6371000;
    p1 = deg2rad(lat1);  p2 = deg2rad(lat2);  dl = deg2rad(lon2 - lon1);
    a = sin((p2 - p1)/2).^2 + cos(p1).*cos(p2).*sin(dl/2).^2;
    d_m = 2*R*atan2(sqrt(a), sqrt(1 - a));
    cog = mod(atan2d(sin(dl).*cos(p2), cos(p1).*sin(p2) - sin(p1).*cos(p2).*cos(dl)), 360);
end

function [yaw, conf] = adaptive_yaw_filter(COG, dist, ok, alpha, refDist)
% Causal circular smoothing of COG; update weight alpha * min(dist/refDist, 1).
% Holds the last estimate when there is no valid observation.
    N = numel(COG);  yaw = nan(N, 1);  conf = zeros(N, 1);
    x = NaN;  y = NaN;  init = false;
    for i = 1:N
        if ok(i) && ~isnan(COG(i))
            w = min(max(dist(i)/refDist, 0), 1);  conf(i) = w;
            if ~init, x = cosd(COG(i));  y = sind(COG(i));  init = true;
            else
                u = alpha*w;
                x = (1 - u)*x + u*cosd(COG(i));  y = (1 - u)*y + u*sind(COG(i));
            end
        end
        if init, yaw(i) = mod(atan2d(y, x), 360); end
    end
end

function y = wrap180(x), y = mod(x + 180, 360) - 180; end
