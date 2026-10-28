%% ocean_wave_spectrum.m
% ground truth for ONE synthetic packet (Holthuijsen description)
% ERA5 sea state at the buoy area 
% synth_truth (PM spectrum, Rayleigh amplitudes, random phases, cos^2s spreading, ice
% attenuation). 

% Figures (same form as the deployment figures):
%   synthetic_target : (a) target E(f) and the realisation, (b) surface elevation
%    synthetic_ice    : (a) open water vs in-ice E(f), (b) energy transmission

% Changes from the previous version:
%   - wdw is used as given (GRIB unit dimensionless = radians); it was passed through
%     deg2rad, which made s hit its cap of 100 (a far too narrow sea).
%   - one fixed random seed (RNG_SEED) so the report figures are reproducible.
%   - the model itself is in synth_truth.m (shared with the 4-day run).

if ~exist('RNG_SEED', 'var'), RNG_SEED = 1; end
rng(RNG_SEED);
FT = sharc_figtools();

%% SETTINGS
ERA5_WAVE = 'Era5Wave.nc';  ERA5_ICE = 'Era5WindIce.nc';
target_lat = -58.8;  target_lon = -0.69;      % Ardhuin coordinates
D        = 1800;     % s, record length (one packet)
CI_FIXED = 1;        % sea-ice concentration used for the attenuation (NaN: use ERA5 siconc)
WCFG = struct('D', D, 'f_max', 0.6, 'ice_dist', 5000, 'a_coeff', 2.12e-3, 'b_coeff', 4.59e-2);
% f_max 0.6 Hz: components above the pipeline's 0.5 Hz anti-alias pass band add nothing
ice_thick = 0.77;    % m, Ardhuin standard thickness (not used by the attenuation model)

%%  ERA5 sea state 
lat  = ncread(ERA5_WAVE, 'latitude');  lon = ncread(ERA5_WAVE, 'longitude');
time = ncread(ERA5_WAVE, 'valid_time');
swh  = ncread(ERA5_WAVE, 'swh');  pp1d = ncread(ERA5_WAVE, 'pp1d');
mwd  = ncread(ERA5_WAVE, 'mwd');  wdw  = ncread(ERA5_WAVE, 'wdw');
[~, ilat] = min(abs(lat - target_lat));  [~, ilon] = min(abs(lon - target_lon));
[~, itime] = max(squeeze(swh(ilon, ilat, :)));        % peak storm in the window
target_time = time(itime);
Hs = swh(ilon, ilat, itime);  Tp = pp1d(ilon, ilat, itime);
mwd_deg = mwd(ilon, ilat, itime);  wdw_v = wdw(ilon, ilat, itime);
ci = CI_FIXED;
if isnan(ci)
    sic = ncread(ERA5_ICE, 'siconc');  ci = sic(ilon, ilat, itime);  if isnan(ci), ci = 0; end
end
fprintf('ERA5 %s: Hs %.2f m, Tp %.1f s, mwd %.0f deg, wdw %.2f rad, ci %.2f\n', ...
    string(datetime(1970,1,1) + seconds(double(target_time))), Hs, Tp, mwd_deg, wdw_v, ci);

%% GROUND TRUTH
W = synth_truth(Hs, Tp, mwd_deg, wdw_v, ci, WCFG);
% unpacked for sharc_meas / direction_validation 
f = W.f;  a = W.a;  a_ice = W.a_ice;  alpha = W.alpha;  theta_i = W.theta_i;
theta_mean = W.theta_mean;  sigma = W.sigma;  s = W.s;  E_f = W.E_f;  Ef_ice = W.Ef_ice;
theta_grid = W.theta_grid;  D_theta = W.D_theta;  transmission = W.transmission;
fp = 1/Tp;
fprintf('spreading: sigma %.0f deg, s = %.1f\n', rad2deg(sigma), s);

% realisation at 2 Hz 
fs = 2;  t = 0:1/fs:D - 1/fs;
eta = zeros(size(t));  eta_ice = zeros(size(t));
for i = 1:numel(f)
    eta     = eta     + a(i)     * cos(2*pi*f(i)*t + alpha(i));
    eta_ice = eta_ice + a_ice(i) * cos(2*pi*f(i)*t + alpha(i));
end

%% ICE DIAGNOSTICS
Hs_open = 4*sqrt(trapz(f, E_f));  Hs_ice = 4*sqrt(trapz(f, Ef_ice));
energy_loss_pct = 100*(1 - trapz(f, Ef_ice)/trapz(f, E_f));
f_half = NaN;
if max(transmission) - min(transmission) > 1e-6
    f_half = interp1(transmission, f, 0.5, 'linear', 'extrap');
end
fprintf('ci %.2f, ice_dist %.0f m: Hs %.2f m -> %.2f m in ice (%.1f%% energy lost), T_E = 0.5 at %.3f Hz\n', ...
    ci, WCFG.ice_dist, Hs_open, Hs_ice, energy_loss_pct, f_half);
m0 = trapz(f, E_f);  m1 = trapz(f, f.*E_f);  m2 = trapz(f, f.^2.*E_f);
fprintf('target: Hm0 %.2f m, Tm01 %.2f s, Tm02 %.2f s, nu %.2f | realisation 4 std(eta) = %.2f m\n', ...
    4*sqrt(m0), m0/m1, sqrt(m0/m2), sqrt(m0*m2/m1^2 - 1), 4*std(eta));

%% TARGET SPCTRUM REALISATION
[Pxx, f_est] = pwelch(eta, hamming(1024), 512, [], fs);
fig = FT.mkfig([900 600]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
a1 = nexttile(tl);
semilogy(a1, f_est, Pxx, '-', 'Color', FT.cb, 'LineWidth', FT.FIG.lw, 'DisplayName', 'Realisation (Welch)');
hold(a1, 'on');  grid(a1, 'on');  box(a1, 'on');
semilogy(a1, f, E_f, '-', 'Color', [0 0 0], 'LineWidth', 0.75, 'DisplayName', 'Target E(f)');
yl = ylim(a1);  ylim(a1, [max(yl(1), 1e-6*max(E_f)) yl(2)]);
xlim(a1, [0 0.5]);  xlabel(a1, 'Frequency (Hz)');  ylabel(a1, 'Spectral Density (m^2/Hz)');
legend(a1, 'Location', 'northeast');
title(a1, sprintf('Target (H_s %.2f m, T_p %.1f s, From %.0f^\\circ): (a) Spectrum, (b) Surface Elevation', ...
    Hs, Tp, mwd_deg), 'FontWeight', 'normal');
a2 = nexttile(tl);  hold(a2, 'on');  grid(a2, 'on');  box(a2, 'on');
iw_ = t <= 300;
plot(a2, t(iw_), eta(iw_), '-', 'Color', FT.cb, 'LineWidth', FT.FIG.lw);
xlabel(a2, 'Time (s)');  ylabel(a2, 'Surface Elevation \eta (m)');
FT.save(fig, 'synthetic_target');

%% ICE ATTENUATION
fig = FT.mkfig([900 600]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
a1 = nexttile(tl);
semilogy(a1, f, E_f, '-', 'Color', [0 0 0], 'LineWidth', 0.75, 'DisplayName', 'Open Water');
hold(a1, 'on');  grid(a1, 'on');  box(a1, 'on');
semilogy(a1, f, Ef_ice, '-', 'Color', FT.cb, 'LineWidth', FT.FIG.lw, ...
    'DisplayName', sprintf('After %.0f km of Ice', WCFG.ice_dist/1e3));
ylim(a1, [1e-6*max(E_f) 10^ceil(log10(max(E_f)))]);
ylabel(a1, 'Spectral Density (m^2/Hz)');  legend(a1, 'Location', 'northeast');  a1.XTickLabel = [];
title(a1, sprintf('Target Spectrum: (a) Open Water and In Ice (c_i = %.2f), (b) Energy Transmission', ci), ...
    'FontWeight', 'normal');
a2 = nexttile(tl);  hold(a2, 'on');  grid(a2, 'on');  box(a2, 'on');
plot(a2, f, transmission, '-', 'Color', FT.cb, 'LineWidth', FT.FIG.lw);
yline(a2, 0.5, ':', 'Color', FT.cg);
if isfinite(f_half), xline(a2, f_half, ':', sprintf('T_E = 0.5 at %.3f Hz', f_half), 'Color', FT.cg, 'FontSize', FT.FIG.fs - 2); end
ylim(a2, [0 1.05]);  ylabel(a2, 'Energy Transmission T_E(f) (-)');  xlabel(a2, 'Frequency (Hz)');
linkaxes([a1 a2], 'x');  xlim(a1, [0 0.5]);
FT.save(fig, 'synthetic_ice');
