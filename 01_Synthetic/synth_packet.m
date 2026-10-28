function [Tpk, S] = synth_packet(W, regime, t0, cfg)
% Raw IMU packet (counts, field CSV columns: DateTimeUTC, Surge, Sway, Heave, Roll, Pitch, Yaw)

%   surface  eta = sum a cos(w t + alpha); slope s = a k khat sin(w t + alpha), khat = travel
%            direction, k = w^2/g (open water dispersion)
%   gyro     pitch rate = d(s_N)/dt, roll rate = d(s_E)/dt (deg/s), yaw 0; + Allan white
%            noise + turn-on bias

%   accel    specific force F = (xdd_E, xdd_N, g + etadd) rotated into the tilted body; the
%            orbital acceleration xdd = -g s cancels g*tilt to first order, so the
%            accelerometer sees almost no tilt at wave frequencies (slopes come from the gyro)
%   counts   accel / scale_heave, gyro / gyro_scale, rounded

% cfg: .fs .T_rec .scale_heave .gyro_scale .n_acc_ug .n_gyro_dps .gyro_bias [p r y]
% S: t, eta, sE, sN (truth signals, for plots)
    g = 9.81;
    amp = W.a(:);  if regime == "ice", amp = W.a_ice(:); end
    fr = W.f(:);  ph0 = W.alpha(:);
    wv = 2*pi*fr;  kv = wv.^2 / g;
    tp = W.theta_i(:) + pi;  cE = cos(tp);  cN = sin(tp);
    t = 0 : 1/cfg.fs : cfg.T_rec - 1/cfg.fs;  nT = numel(t);
    [eta, etadd, sE, sN, dsE, dsN, xddE, xddN] = deal(zeros(1, nT));
    CH = 25;
    for i0 = 1:CH:numel(fr)
        j = i0:min(i0 + CH - 1, numel(fr));
        P = wv(j) * t + ph0(j);  C = cos(P);  Sn = sin(P);
        A = amp(j);  Wj = wv(j);  K = kv(j);
        eta   = eta   + sum(A .* C, 1);
        etadd = etadd - sum(A .* Wj.^2 .* C, 1);
        sE    = sE    + sum(A .* K .* cE(j) .* Sn, 1);
        sN    = sN    + sum(A .* K .* cN(j) .* Sn, 1);
        dsE   = dsE   + sum(A .* K .* Wj .* cE(j) .* C, 1);
        dsN   = dsN   + sum(A .* K .* Wj .* cN(j) .* C, 1);
        xddE  = xddE  - sum(A .* Wj.^2 .* cE(j) .* Sn, 1);
        xddN  = xddN  - sum(A .* Wj.^2 .* cN(j) .* Sn, 1);
    end
    th = sN;  ph = sE;  Fup = g + etadd;
    fx = xddN .* cos(th) + Fup .* sin(th);
    fy = xddE .* cos(ph) + Fup .* sin(ph);
    fz = -xddE .* sin(ph) - xddN .* sin(th) + Fup .* cos(th) .* cos(ph);
    sa = cfg.n_acc_ug*1e-6*9.80665 * sqrt(cfg.fs/2);
    sg = cfg.n_gyro_dps * sqrt(cfg.fs/2);
    n_ = @() randn(1, nT);
    Surge = round((fx + sa*n_()) / cfg.scale_heave);
    Sway  = round((fy + sa*n_()) / cfg.scale_heave);
    Heave = round((fz + sa*n_()) / cfg.scale_heave);
    Pitch = round((rad2deg(dsN) + sg*n_() + cfg.gyro_bias(1)) / cfg.gyro_scale);
    Roll  = round((rad2deg(dsE) + sg*n_() + cfg.gyro_bias(2)) / cfg.gyro_scale);
    Yaw   = round((sg*n_() + cfg.gyro_bias(3)) / cfg.gyro_scale);
    DateTimeUTC = t0 + seconds(t(:));
    Tpk = table(DateTimeUTC, Surge(:), Sway(:), Heave(:), Roll(:), Pitch(:), Yaw(:), ...
        'VariableNames', {'DateTimeUTC','Surge','Sway','Heave','Roll','Pitch','Yaw'});
    S = struct('t', t, 'eta', eta, 'sE', sE, 'sN', sN);
end
