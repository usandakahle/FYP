%% sharc_meas.m -- ONE synthetic packet through the SAME pipeline as the field data
% Raw IMU packet (field CSV columns, counts) from the ground truth in ocean_wave_spectrum.m
% (synth_packet.m), processed by run_sharc_jack with the options batchBias uses. No alpha_c:
% tilt fusion is the pipeline's Kalman/PI filter (kP, kI from the Allan variance).
% Packet figure S1 in exactly the field form (shared code: sharc_figtools packet_fig).
% For the 4-day run (20-24 July) see synthetic_deployment.m.

if ~exist('W', 'var')
    run('ocean_wave_spectrum.m');          % ground truth W (and the unpacked variables)
end

%% ---------------- CONFIG ----------------
if ~exist('REGIME', 'var'), REGIME = "ice"; end      % "open" | "ice"
REGIME = string(REGIME);
PCFG = struct('fs', 100, 'T_rec', D, 'scale_heave', 9.80665/8192, 'gyro_scale', 1/65.5, ...
    'n_acc_ug', 593, 'n_gyro_dps', 0.0352, 'gyro_bias', [0.3 -0.2 0.1]);
% noise = bench Allan white noise (same numbers as the pipeline noise model); gyro bias is
% removed by the pipeline (burst median)
opts = struct('fs_imu','auto', 'gyro_scale_dps',PCFG.gyro_scale, 'n_acc_ug',PCFG.n_acc_ug, ...
    'n_gyro_dps',PCFG.n_gyro_dps, 'snr_min',3, 'noise_correct',true, 'f_min',0.03, 'f_max_band',0.35, ...
    'Nfft',512, 'M',32, 'tau',37, 'theta_offset_deg',0, 'yaw_sign',1, 'gnss_yaw_sense',1, ...
    'heading_offset_deg',0, 'jackknife',false, 'track',false, 'verbose',true, 'gnss_yaw',[]);
% same opts as batchBias; gnss_yaw = [] because the synthetic bow stays north (body = earth)
SAVE_PACKET = false;                 % true: also write the packet CSV in the field format
POLAR_FMAX  = 0.30;
FT = sharc_figtools();

%% ---------------- packet -> pipeline ----------------
t0 = datetime(1970,1,1, 'TimeZone','UTC') + seconds(double(target_time));
[Tpk, SIG] = synth_packet(W, REGIME, t0, PCFG);
if SAVE_PACKET
    if ~exist(FT.FIG.outdir, 'dir'), mkdir(FT.FIG.outdir); end
    writetable(Tpk, fullfile(FT.FIG.outdir, "packet_synthetic_V2" + REGIME + ".csv"));
end
[theta_sharc, Hm0_sharc, Tm01_sharc, ~, pd] = run_sharc_jack(Tpk, "synthetic " + REGIME, false, ...
    PCFG.scale_heave, 0, +1, [], opts);
R = synth_compare(W, REGIME, pd, opts.f_min);

% names used by direction_validation.m
dir_target = R.dir_target;  dir_real = R.dir_real;  Czc = R.Czc;  chk_f = R.chk_f;
fr = R.f;  E_tgt = R.E_target;  dfp = pd.freq(2) - pd.freq(1);
t_imu = SIG.t;  eta_t = SIG.eta;  sE = SIG.sE;  sN = SIG.sN;

fprintf('\n=== synthetic packet (%s), same pipeline as the field data ===\n', REGIME);
fprintf('%-22s %9s %9s %9s\n', '', 'target', 'realised', 'SHARC');
fprintf('%-22s %9.2f %9.2f %9.2f\n', 'Hm0 (m)', R.Hm0_target, R.Hm0_real, Hm0_sharc);
fprintf('%-22s %9.2f %9.2f %9.2f\n', 'Tm01 (s)', R.Tm01_target, R.Tm01_real, Tm01_sharc);
fprintf('%-22s %9.1f %9.1f %9.1f   (error vs realised %+.1f)\n', 'direction FROM (deg)', ...
    dir_target, dir_real, theta_sharc, FT.wrap180(theta_sharc - dir_real));
fprintf('%-22s %9.1f %9s %9.1f\n', 'spread (deg)', rad2deg(sigma), '', pd.spread_sys);
fprintf('coherence %.2f, check ratio %.2f, valid %d %s\n', pd.coh_sys, R.chk_sys, pd.dir_valid, pd.invalid_reason);

%% ---------------- S1: packet figure (field form) ----------------
rlab = "Open Water";  if REGIME == "ice", rlab = "In Ice"; end
head = {sprintf('Synthetic Packet (%s)', rlab), ...
    sprintf('H_{m0} = %.2f m,   \\theta_{alg} = %.0f^\\circ,   Target = %.0f^\\circ,   r_1 = %.2f', ...
    Hm0_sharc, theta_sharc, dir_target, pd.r1_sys)};
dirs = struct('deg', {theta_sharc, dir_target}, 'col', {FT.cb, [0 0 0]}, 'ls', {'-', '--'}, ...
    'name', {'SHARC', 'Target'});
fig = FT.packet_fig(pd, head, dirs, POLAR_FMAX, struct('f', fr, 'E', E_tgt, 'name', 'Target E(f)'));
FT.save(fig, char("synthetic_packet_" + REGIME));

%% ---------------- attenuation estimate from the SHARC spectrum (unchanged method) ----------------
fa_ = pd.freq(:);  low = fa_ > opts.f_min & fa_ < 0.3 & Czc > 0;
pm  = @(p, ff) 0.3125*p(1)^2*(1/p(2))^4 ./ ff.^5 .* exp(-1.25*(1./(p(2)*ff)).^4);   % p = [Hs Tp]
[~, ipk] = max(Czc .* low);
p_hat = fminsearch(@(p) sum((pm(p, fa_(low)) - Czc(low)).^2), [Hm0_sharc, 1/fa_(ipk)]);
fit_b = fa_ > 1/p_hat(2) & fa_ < 0.5 & Czc > 0;
AB = [-fa_(fit_b).^2, -fa_(fit_b).^4] \ log(max(Czc(fit_b) ./ pm(p_hat, fa_(fit_b)), eps));
fprintf('PM fit: Hs %.2f m, Tp %.1f s | attenuation fit A = %.3g s^2, B = %.3g s^4\n', p_hat, AB);
