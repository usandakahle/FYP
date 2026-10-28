function W = synth_truth(Hs, Tp, mwd_deg, wdw, ci, cfg)
% Ground-truth wave field for one packet (Holthuijsen random-phase model, as in
% ocean_wave_spectrum.m). Used by ocean_wave_spectrum.m (one packet) and
% synthetic_deployment.m (20-24 July), so both use the same model.
%   Hs (m), Tp (s), mwd_deg: ERA5 mean direction FROM (deg, clockwise from N)
%   wdw: ERA5 wave spectral directional width. GRIB unit is dimensionless (a directional
%        standard deviation in radians, ~0.3-1.4); values > pi are taken as degrees.
%   ci : sea-ice concentration (0-1), scales the attenuation (as before)
%   cfg: .D (s) record length, .f_max (Hz) highest component, .ice_dist (m),
%        .a_coeff, .b_coeff (Meylan 2014 / Doble 2015 pancake coefficients)
% W: f, a, a_ice, alpha, theta_i (math angle of FROM), theta_mean, sigma (rad), s, Cs,
%    E_f, Ef_ice, transmission, Hs, Tp, mwd_deg, ci
    Hs = double(Hs);  Tp = double(Tp);  mwd_deg = double(mwd_deg);  wdw = double(wdw);  ci = double(ci);
    W.Hs = Hs;  W.Tp = Tp;  W.mwd_deg = mwd_deg;  W.ci = ci;          % (ncread gives single)
    df = 1/cfg.D;
    f  = df:df:cfg.f_max;  N_f = numel(f);
    fp = 1/Tp;
    E_f = (0.3125 * Hs^2 * fp^4 ./ f.^5) .* exp(-1.25 * (fp./f).^4);    % Pierson-Moskowitz
    alpha = 2*pi*rand(1, N_f);
    a = raylrnd(sqrt(E_f*df));                                        % Rayleigh amplitudes
    % attenuation in ice (Meylan 2014 / Doble 2015), scaled by ci
    alpha_eff = ci * (cfg.a_coeff*f.^2 + cfg.b_coeff*f.^4);
    W.transmission = exp(-alpha_eff * cfg.ice_dist);
    W.Ef_ice = E_f .* W.transmission;
    W.a_ice  = a .* exp(-0.5*alpha_eff*cfg.ice_dist);
    % directional spreading cos^(2s), one direction per component
    sigma = wdw;  if sigma > pi, sigma = deg2rad(sigma); end
    s = min(max(2/sigma^2 - 1, 1), 100);
    logCs = (2*s-1)*log(2) - log(pi) + 2*gammaln(s+1) - gammaln(2*s+1);
    theta_mean = wrapToPi(deg2rad(90 - mwd_deg));                    % math angle of FROM
    tg = linspace(-pi, pi, 721);
    Dth = exp(logCs) * max(cos(wrapToPi(tg - theta_mean)/2), 0).^(2*s);
    cdf = cumtrapz(tg, Dth);  cdf = cdf / cdf(end);
    [cu, iu] = unique(cdf);
    W.theta_i = interp1(cu, tg(iu), rand(1, N_f), 'linear', 'extrap');
    W.f = f;  W.a = a;  W.alpha = alpha;  W.E_f = E_f;
    W.theta_mean = theta_mean;  W.sigma = sigma;  W.s = s;  W.Cs = exp(logCs);
    W.theta_grid = tg;  W.D_theta = Dth;
end
