%% SHARC buoy measurement
% Simulates the SHARC sensing and DSP chain from a known wave spectrum.
% Ground truth is generated in ocean_wave_spectrum.m.

clear; clc;
run('ocean_wave_spectrum.m');

%% 1. Model settings
% D       : record length; controls the number of Bartlett segments.
%           A short record gives a noisy PSD and unreliable Hm0/Tm01/direction.
% f_min   : removes low-frequency drift before acceleration is converted
%           to displacement. Double integration amplifies low frequencies as 1/f^2.
% fs_imu  : IMU sample rate; must be high enough for the wave band.
% noise_density / gyro_noise_density : sensor noise floors.
% alpha_c : complementary-filter weighting between gyro and accelerometer.
% M       : decimation factor; determines the final sample rate.
% cutoff_norm : FIR low-pass cutoff; too low removes real wave energy.
% tau     : RC high-pass time constant; removes drift but can attenuate swell.
% Nfft    : FFT segment length; sets frequency resolution fs_dec/Nfft.
% window  : controls spectral leakage versus frequency resolution.


%% Heave Acceleration
% eta(t) = sum_i a(i)*cos(2*pi*f(i)*t + alpha(i))
% Differentiate the same sum analytically (avoids numerical-diff noise):
% d^2/dt^2 [a*cos(wt+phi)] = -a*w^2*cos(wt+phi)
%======================================================
g = 9.81;
w = 2*pi*f;                      
k = w.^2 / g;                     %Deep water wave number
fs_imu = 100;                     % (Hz) ICM20649 sample rate
t_imu = 0: 1/fs_imu: D-1/fs_imu;   %time vector for imu
dt_imu = 1/fs_imu;

%% Eta
% vertical accelartion as the second derivative of the surface elevation
%η(t) = Σ a·cos(wt+φ)
%===============================================================
heave_spec = zeros(size(t_imu));
for i = 1:length(f)
    heave_spec = heave_spec - a(i)*w(i)^2*cos(w(i)*t_imu + alpha(i));
end


%% Roll/Pitch
% components of the slope (rate of change of the slope) (deriv)
%Sx = -a·k·cos(θ)·sin(wt+φ)
%Sy = -a·k·sin(θ)·sin(wt+φ)
% theta os the direction of the wave components (much like a pythagorean
% triangle)
%=====================================================
gx = zeros(size(t_imu));   % true pitch rate (rad/s)
gy = zeros(size(t_imu));   % true roll rate (rad/s)
for i = 1:length(f)
    gx = gx - a(i)*k(i)*cos(theta_i(i))*w(i)*cos(w(i)*t_imu + alpha(i));
    gy = gy - a(i)*k(i)*sin(theta_i(i))*w(i)*cos(w(i)*t_imu + alpha(i));
end

%% True slope (for reference)
% The slopes that we attempt to estimate (for comparison of how accurate
% the model is
%===================================================
Sx_true = zeros(size(t_imu));
Sy_true = zeros(size(t_imu));
for i = 1:length(f)
    Sx_true = Sx_true - a(i)*k(i)*cos(theta_i(i)) * sin(w(i)*t_imu + alpha(i));
    Sy_true = Sy_true - a(i)*k(i)*sin(theta_i(i)) * sin(w(i)*t_imu + alpha(i));
end

%% Gyro model
%======================================================
gyro_noise_density = 0.0175 * (pi/180);          % deg/s/sqrt(Hz) -> rad/s/sqrt(Hz)
sigma_gyro = gyro_noise_density * sqrt(fs_imu/2);
gx_meas = gx + sigma_gyro * randn(size(t_imu));
gy_meas = gy + sigma_gyro * randn(size(t_imu));

%% IMU sensor model (ICM20649)
%======================================================================

noise_density = 289e-6 * 9.81;   % 289 ug/sqrt(Hz) -> m/s^2/sqrt(Hz)
sigma_noise = noise_density * sqrt(fs_imu/2);
accel_noise = sigma_noise * randn(size(t_imu));
% 32-bit stored, 4g full scale: quantisation step is negligible at these
% amplitudes (waves-in-ice < 1 m/s^2), so it's omitted here.

heave = heave_spec + accel_noise;


%% Accelerometer-derived static tilt (small-angle: horizontal accel ~= g*slope)
%small angle aprox, simplifies
%Angular tilt of the IMU/buoy
%We assume the buoy is rigidly attaced to the stand, and doesnt rotate such
%that the tilt of the stand is approx equal to the tilt of the buoy
%===================================================================================
ax_tilt_meas = g*Sx_true + sigma_noise * randn(size(t_imu));  % sigma_noise from heave IMU block
ay_tilt_meas = g*Sy_true + sigma_noise * randn(size(t_imu));
Sx_accel_est = ax_tilt_meas / g;
Sy_accel_est = ay_tilt_meas / g;

%% Complementary filter: gyro (short-term) + accel tilt (long-term drift correction)
%gyro correction
%=================================================================================
alpha_c = 0.98;    % blend constant -- higher trusts the gyro more, corrects slowly from accel
%design choice
Sx_fused = zeros(size(t_imu));
Sy_fused = zeros(size(t_imu));
for n = 2:length(t_imu)
    Sx_fused(n) = alpha_c*(Sx_fused(n-1) + gx_meas(n)*dt_imu) + (1-alpha_c)*Sx_accel_est(n);
    Sy_fused(n) = alpha_c*(Sy_fused(n-1) + gy_meas(n)*dt_imu) + (1-alpha_c)*Sy_accel_est(n);
end


%% Already Existing DSP pipeline 
% Detrend raw block(removes mean data)
%===============================================================
heave= detrend(heave);                  %removes linear trend

% FIR low-pass + decimate, M = 32 (100 Hz -> 3.125 Hz nominal)
M = 32;                                   %change
fs_dec = fs_imu / M;
if fs_imu/2 <= fs_dec
    warning(['fs_imu = %.2f Hz cannot be meaningfully decimated by M=%d ' ...
        '(Nyquist is already at/under the target rate). Numbers below ' ...
        'are structural only — see the mismatch note at the top.'], fs_imu, M);
end

%% Anti-Aliasing
%Drived from thesis
%==================================================
cutoff_norm = min(0.9, 1.12 / (fs_imu/2));   % thesis cutoff = 1.12 Hz
lpf = fir1(64, cutoff_norm);                 %65 taps
%change and observe changes sharpness vs latency
heave_filt = filtfilt(lpf, 1, heave);          %reconsider (not causal so may not work in real time)
heave_dec = downsample(heave_filt, M);

%  RC high-pass detrend, tau = 100 s (thesis: k=0.9986 @ 6.25 Hz;
% same formula generalised to whatever fs_dec is found t \o be
tau = 100;

k_rc = exp(-1/(tau*fs_dec));

heave_hp = zeros(size(heave_dec));
for n = 2:length(heave_dec)
    heave_hp(n) = k_rc*heave_hp(n-1) + k_rc*(heave_dec(n) - heave_dec(n-1));
end

%% FFT Segmentation
%========================================================
% Segment into 1024-sample blocks for FFT (Cortex-M4 CMSIS cap)
Nfft = 1024;
Nsamples = length(heave_hp);
Nseg = floor(Nsamples / Nfft);
if Nseg < 1
    error(['Only %d decimated samples available — need at least %d for one ' ...
        'FFT segment. Increase D in ocean_wave_spectrum.m.'], Nsamples, Nfft);
elseif Nseg < 6
    warning(['Only %d full %d-sample segment(s) available; the thesis ' ...
        'Bartlett average uses N=6. Increase record length for a closer match.'], ...
        Nseg, Nfft);
end

%% Windowing Function
%======================================================
win = hanning(Nfft);                   %windowing function
%win_energy_correction = 1.63;      % thesis-quoted Hanning energy correction
freq_axis = (0:Nfft/2-1) * (fs_dec/Nfft);
%freq_axis(1) = freq_axis(2)/2;     % avoid divide-by-zero at DC below

%% Segments in the PSD (one sided)
seg_psds = zeros(Nfft/2, Nseg);
for s = 1:Nseg                 % for each segment
    idx = (s-1)*Nfft + (1:Nfft);
    x = heave_hp(idx) .* win(:)';

    X = fft(x, Nfft);
    X = X(1:Nfft/2);                       

    f_min = 0.03;   % Hz, thesis drift cutoff (~33 s period)
    % acceleration -> displacement: divide by -(2*pi*f)^2
    valid = freq_axis >= f_min & freq_axis <= 0.5;   %upper band for useful plots
    %H = -1 ./ (2*pi*freq_axis).^2;
    X_disp = zeros(size(X));
    X_disp(valid) = X(valid) .* (-1 ./ (2*pi*freq_axis(valid)).^2);

    %  low-pass the displacement spectrum to cut leakage above the
    % wave band (thesis doesn't give an exact cutoff for this stage —
    % 0.5 Hz used here as a placeholder; tighten if you find the real one)
    %X_disp = X_disp .* (freq_axis <= 0.5);

    % 3f. band-limit the displacement spectrum: cut leakage above the wave
    % band AND exclude near-DC bins below the thesis drift cutoff, since
    % those get blown up by 1/f^2 in the double-integration step above.
   
    X_disp = X_disp .* (freq_axis >= f_min & freq_axis <= 0.5);

    %  PSD, energy-corrected, Parseval-scaled
    Pxx = (2 / (fs_dec * sum(win.^2))) * abs(X_disp).^2;
    
    seg_psds(:,s) = Pxx;
end

%% Slope channel DSP (parallel to heave -- same decimate/detrend/window/FFT chain)
%=====================================================================
Sx_filt = filtfilt(lpf, 1, Sx_fused);
Sy_filt = filtfilt(lpf, 1, Sy_fused);
Sx_dec = downsample(Sx_filt, M);
Sy_dec = downsample(Sy_filt, M);

Sx_hp = zeros(size(Sx_dec));
Sy_hp = zeros(size(Sy_dec));
for n = 2:length(Sx_dec)
    Sx_hp(n) = k_rc*Sx_hp(n-1) + k_rc*(Sx_dec(n) - Sx_dec(n-1));
    Sy_hp(n) = k_rc*Sy_hp(n-1) + k_rc*(Sy_dec(n) - Sy_dec(n-1));
end

% NOTE: unlike heave, slope channels are NOT double-integrated
% (Sx_fused/Sy_fused already estimate slope directly, not acceleration),
% so there is no -1/(2*pi*f)^2 factor applied here.

seg_Z = zeros(Nfft/2, Nseg);   % complex heave FFT per segment (for cross-spectra)
seg_X = zeros(Nfft/2, Nseg);   % complex Sx FFT per segment
seg_Y = zeros(Nfft/2, Nseg);   % complex Sy FFT per segment


%% Welch-averaged version (For comparison)
%2 times more segments for smoother PSD
%===========================================================================
overlap_frac = 0.5;                        % 50% overlap = classic Welch
hop = round(Nfft * (1 - overlap_frac));
Nseg_w = floor((Nsamples - Nfft) / hop) + 1;

seg_psds_w = zeros(Nfft/2, Nseg_w);
seg_Z_w = zeros(Nfft/2, Nseg_w);
seg_X_w = zeros(Nfft/2, Nseg_w);
seg_Y_w = zeros(Nfft/2, Nseg_w);

for s = 1:Nseg_w
    idx = (s-1)*hop + (1:Nfft);

    xz = heave_hp(idx) .* win(:)';
    Xz = fft(xz, Nfft); Xz = Xz(1:Nfft/2);
    Xz_disp = zeros(1, Nfft/2);
    Xz_disp(valid) = Xz(valid) .* (-1 ./ (2*pi*freq_axis(valid)).^2);
    seg_psds_w(:,s) = (2 / (fs_dec * sum(win.^2))) * abs(Xz_disp).^2;
    seg_Z_w(:,s) = Xz_disp;

    xx = Sx_hp(idx) .* win(:)';
    Xx = fft(xx, Nfft); Xx = Xx(1:Nfft/2); Xx(~valid) = 0;
    seg_X_w(:,s) = Xx;

    xy = Sy_hp(idx) .* win(:)';
    Xy = fft(xy, Nfft); Xy = Xy(1:Nfft/2); Xy(~valid) = 0;
    seg_Y_w(:,s) = Xy;
end
%% Welch PSD
%===============================
Pxx_avg_w = mean(seg_psds_w, 2);

valid = freq_axis >= f_min & freq_axis <= 0.5;

for s = 1:Nseg
    idx = (s-1)*Nfft + (1:Nfft);

    % --- Heave (re-derive complex spectrum here since the existing loop
    %     only kept abs(X_disp).^2 in seg_psds, not the complex value) ---
    xz = heave_hp(idx) .* win(:)';
    Xz = fft(xz, Nfft); Xz = Xz(1:Nfft/2);
    Xz_disp = zeros(1, Nfft/2);
    Xz_disp(valid) = Xz(valid) .* (-1 ./ (2*pi*freq_axis(valid)).^2);
    seg_Z(:,s) = Xz_disp;

    % --- Sx, Sy (no acceleration->displacement inversion needed) ---
    xx = Sx_hp(idx) .* win(:)';
    Xx = fft(xx, Nfft); Xx = Xx(1:Nfft/2);
    Xx(~valid) = 0;
    seg_X(:,s) = Xx;

    xy = Sy_hp(idx) .* win(:)';
    Xy = fft(xy, Nfft); Xy = Xy(1:Nfft/2);
    Xy(~valid) = 0;
    seg_Y(:,s) = Xy;
end

%%  Bartlett-averaged PSD (consider welch)
%==========================================================
Pxx_avg = mean(seg_psds, 2);


%%  Wave parameters from the SHARC-measured PSD 
%=======================================================================
m0 = trapz(freq_axis, Pxx_avg);                     %variance of surface elevation
m1 = trapz(freq_axis, freq_axis .* Pxx_avg);         %first moment
Hm0_sharc = 4*sqrt(m0);                             %significant wave height    
Tm01_sharc = m0./m1                                 %mean period (not peak)

fprintf('SHARC-measured : Hm0 = %.3f m, Tm01 = %.2f s\n', Hm0_sharc, Tm01_sharc);
fprintf('Input truth    : Hs  = %.3f m, Tp  = %.2f s (peak, not Tm01 — compare with care)\n', Hs, Tp);

%%  Plots
%==========================================================
%{figure;

%plot(t, heave_spec, t, heave, '--');
%xlabel('Time (s)'); ylabel('Heave accel (m/s^2)');
%legend('true heave','IMU-measured heave');
%title('Heave acceleration: truth vs sensor');

%{
h1 = plot(t_imu, heave_spec, 'LineWidth', 1.2);
hold on;
h2 = plot(t_imu, heave, '--', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Heave accel (m/s^2)');
lgd = legend([h1 h2], {'True heave','IMU-measured heave'});
title('Heave acceleration: truth vs sensor');
grid on;

% Make legend entries clickable
lgd.ItemHitFcn = @(~,event) set(event.Peer, 'Visible', ~strcmp(event.Peer.Visible,'on') );
%}

figure;
h3 = plot(f, E_f,'LineWidth',1.2);
hold on;
h4 = plot(freq_axis, Pxx_avg, 'LineWidth', 1.2);
xlabel('Frequency (Hz)'); 
ylabel('Displacement PSD (m^2/Hz)');
lgd = legend([h3 h4], {'Target Spectrum','IMU-measured Spectrum'});
title('SHARC-recovered heave displacement spectrum');
grid on;

%{
%% Plot A: Gyro measurement vs truth (short window)
idx_win = t_imu <= 100;
figure;
subplot(2,1,1);
plot(t_imu(idx_win), gx(idx_win), 'LineWidth', 1.2); hold on;
plot(t_imu(idx_win), gx_meas(idx_win), '--', 'LineWidth', 1.0);
ylabel('Pitch rate (rad/s)'); legend('True g_x', 'Measured g_x'); grid on;
title('Gyro measurement vs truth (first 100 s)');
subplot(2,1,2);
plot(t_imu(idx_win), gy(idx_win), 'LineWidth', 1.2); hold on;
plot(t_imu(idx_win), gy_meas(idx_win), '--', 'LineWidth', 1.0);
xlabel('Time (s)'); ylabel('Roll rate (rad/s)'); legend('True g_y', 'Measured g_y'); grid on;

%% Plot B: Complementary filter fusion vs truth (does drift correction work?)
figure;
subplot(2,1,1);
plot(t_imu(idx_win), Sx_true(idx_win), 'LineWidth', 1.5); hold on;
plot(t_imu(idx_win), Sx_fused(idx_win), 'LineWidth', 1.0);
plot(t_imu(idx_win), Sx_accel_est(idx_win), ':', 'LineWidth', 0.8);
ylabel('S_x (rad)'); legend('True', 'Fused (comp. filter)', 'Accel-only estimate'); grid on;
title('Slope fusion: truth vs complementary filter (first 100 s)');
subplot(2,1,2);
plot(t_imu(idx_win), Sy_true(idx_win), 'LineWidth', 1.5); hold on;
plot(t_imu(idx_win), Sy_fused(idx_win), 'LineWidth', 1.0);
plot(t_imu(idx_win), Sy_accel_est(idx_win), ':', 'LineWidth', 0.8);
xlabel('Time (s)'); ylabel('S_y (rad)'); legend('True', 'Fused', 'Accel-only'); grid on;

% If Sx_fused visibly drifts away from Sx_true over a LONGER window than
% shown here, alpha_c is too high (trusting gyro drift too long) --
% check on a 300+ s window before trusting the DSP-chain results below.
figure;
idx_long = t_imu <= min(600, D);
plot(t_imu(idx_long), Sx_true(idx_long), 'LineWidth', 1.2); hold on;
plot(t_imu(idx_long), Sx_fused(idx_long), 'LineWidth', 1.0);
xlabel('Time (s)'); ylabel('S_x (rad)');
legend('True', 'Fused'); title('Long-window drift check (S_x)'); grid on;
%}

%% Plot C: Fused slope PSD vs theoretical (post-fusion, pre-decimation check)
[Pxx_Sx_fused, f_fused] = pwelch(Sx_fused, hamming(1024), 512, [], fs_imu);
[Pxx_Sy_fused, ~]       = pwelch(Sy_fused, hamming(1024), 512, [], fs_imu);



figure;
subplot(2,1,1);
%plot(f, E_Sx_theory, 'LineWidth', 1.5); hold on;
plot(f_fused, Pxx_Sx_fused, 'LineWidth', 1.0);
xlim([0 1]); ylabel('S_x PSD (rad^2/Hz)');
legend('Theoretical', 'Fused (pre-decimation)'); title('S_x PSD: fusion output vs theory'); grid on;
subplot(2,1,2);
%plot(f, E_Sy_theory, 'LineWidth', 1.5); hold on;
plot(f_fused, Pxx_Sy_fused, 'LineWidth', 1.0);
xlim([0 1]); xlabel('Frequency (Hz)'); ylabel('S_y PSD (rad^2/Hz)');
legend('Theoretical', 'Fused (pre-decimation)'); grid on;


%{
%% Plot D: Decimated/detrended slope channels (mirrors existing heave diagnostic plots)
figure;
subplot(3,1,1);
plot((0:length(Sx_dec)-1)/fs_dec, Sx_dec);
title('Decimated S_x'); ylabel('S_x (rad)'); grid on;
subplot(3,1,2);
plot((0:length(Sx_hp)-1)/fs_dec, Sx_hp);
title('Detrended (RC high-pass) S_x'); ylabel('S_x (rad)'); grid on;
subplot(3,1,3);
plot((0:length(Sy_hp)-1)/fs_dec, Sy_hp);
title('Detrended (RC high-pass) S_y'); xlabel('Time (s)'); ylabel('S_y (rad)'); grid on;

%} 

%% Plot E: Bartlett-averaged auto-spectra of z, Sx, Sy (sanity before cross-spectra)
K_norm = 2 / (fs_dec * sum(win.^2));   % same normalisation as existing Pxx calculation

Czz_avg = K_norm * mean(abs(seg_Z).^2, 2);
Cxx_avg = K_norm * mean(abs(seg_X).^2, 2);
Cyy_avg = K_norm * mean(abs(seg_Y).^2, 2);

figure;
plot(freq_axis, Czz_avg, 'LineWidth', 1.2); hold on;
plot(freq_axis, Cxx_avg, 'LineWidth', 1.2);
plot(freq_axis, Cyy_avg, 'LineWidth', 1.2);
xlim([0 0.3]);
xlabel('Frequency (Hz)'); ylabel('PSD');
legend('C_{zz} (heave)', 'C_{xx} (S_x)', 'C_{yy} (S_y)');
title('Bartlett-averaged auto-spectra, SHARC pipeline output'); grid on;


%% Plot F: Directional recovery -- the actual result
Qzx_avg = -K_norm * mean(imag(seg_Z .* conj(seg_X)), 2);
Qzy_avg = -K_norm * mean(imag(seg_Z .* conj(seg_Y)), 2);
Cxy_avg = K_norm * mean(real(seg_X .* conj(seg_Y)), 2);

%% Directional Fourier Coefficient
a1 = Qzx_avg ./ sqrt(Czz_avg .* (Cxx_avg + Cyy_avg));     
b1 = Qzy_avg ./ sqrt(Czz_avg .* (Cxx_avg + Cyy_avg));
a2 = (Cxx_avg - Cyy_avg) ./ (Cxx_avg + Cyy_avg);
b2 = 2*Cxy_avg ./ (Cxx_avg + Cyy_avg);

%theta1_deg = rad2deg(wrapTo2Pi(atan2(b1, a1)));
theta1_deg = mod(90 - rad2deg(atan2(b1, a1)), 360);   % math -> compass bearing
spread_deg = rad2deg(sqrt(max(2*(1 - sqrt(a1.^2 + b1.^2)), 0)));

band = freq_axis >= f_min & freq_axis <= 0.3;   % restrict to the physically valid band


figure;
pax = polaraxes;
polarplot(pax, deg2rad(theta1_deg(band)), freq_axis(band), 'o', 'MarkerSize', 4); hold on;
%polarplot(pax, theta_mean*ones(1,2), [0 max(freq_axis(band))], 'r--', 'LineWidth', 1.5);

theta_mean_compass_deg = mod(90 - rad2deg(theta_mean), 360);
polarplot(pax, deg2rad(theta_mean_compass_deg)*ones(1,2), [0 max(freq_axis(band))], 'r--', 'LineWidth', 1.5);
pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
title('Recovered direction per frequency bin (radius = frequency)');


figure;
subplot(2,1,1);
plot(freq_axis(band), theta1_deg(band), 'o-', 'LineWidth', 1.0); hold on;
%yline(rad2deg(wrapTo2Pi(theta_mean)), 'r--', 'LineWidth', 1.5);
yline(theta_mean_compass_deg, 'r--', 'LineWidth', 1.5);
ylabel('\theta_1(f) (deg)'); legend('Recovered', 'Input \theta_{mean}'); grid on;
title('SHARC-recovered direction and spread vs frequency');
subplot(2,1,2);
plot(freq_axis(band), spread_deg(band), 'o-', 'LineWidth', 1.0);
xlabel('Frequency (Hz)'); ylabel('Spread (deg)'); grid on;

%% Welch auto/cross-spectra + directional coefficients
Czz_avg_w = K_norm * mean(abs(seg_Z_w).^2, 2);
Cxx_avg_w = K_norm * mean(abs(seg_X_w).^2, 2);
Cyy_avg_w = K_norm * mean(abs(seg_Y_w).^2, 2);

Qzx_avg_w = -K_norm * mean(imag(seg_Z_w .* conj(seg_X_w)), 2);
Qzy_avg_w = -K_norm * mean(imag(seg_Z_w .* conj(seg_Y_w)), 2);
Cxy_avg_w = K_norm * mean(real(seg_X_w .* conj(seg_Y_w)), 2);

a1_w = Qzx_avg_w ./ sqrt(Czz_avg_w .* (Cxx_avg_w + Cyy_avg_w));
b1_w = Qzy_avg_w ./ sqrt(Czz_avg_w .* (Cxx_avg_w + Cyy_avg_w));
a2_w = (Cxx_avg_w - Cyy_avg_w) ./ (Cxx_avg_w + Cyy_avg_w);
b2_w = 2*Cxy_avg_w ./ (Cxx_avg_w + Cyy_avg_w);

%% Reconstruct E(f,theta) for a given Czz, a1,b1,a2,b2
% theta_grid below IS the compass bearing used for both the physics and the plot geometry.
build_E2D_compass = @(Czz,a1,b1,a2,b2,theta_c_rad) ...
    Czz .* max( (1/pi)*(0.5 + a1.*sin(theta_c_rad) + b1.*cos(theta_c_rad) ...
                 - a2.*cos(2*theta_c_rad) + b2.*sin(2*theta_c_rad)), 0 );

theta_grid = (0:2:358)';              % compass bearing, deg -- used for BOTH physics and geometry now
theta_c_rad = deg2rad(theta_grid)';   % row vector

E2D_bart  = build_E2D_compass(Czz_avg,   a1,   b1,   a2,   b2,   theta_c_rad);
E2D_welch = build_E2D_compass(Czz_avg_w, a1_w, b1_w, a2_w, b2_w, theta_c_rad);

band = freq_axis >= f_min & freq_axis <= 0.3;
f_band = freq_axis(band);

E2D_bart_norm  = E2D_bart(band,:)  / max(E2D_bart(band,:), [], 'all');
E2D_welch_norm = E2D_welch(band,:) / max(E2D_welch(band,:), [], 'all');

%% Plot G: Bartlett vs Welch directional spectra — compass-bearing heatmaps
% theta_grid here is already compass bearing (deg)
fmax  = 0.3;
[Theta, F] = meshgrid(theta_grid, f_band);
Xc = F .* sind(Theta);
Yc = F .* cosd(Theta);

figure('Position',[100 100 600 600]);
pcolor(Xc, Yc, E2D_bart_norm); shading interp; colormap(jet); axis equal; hold on;
add_polar_grid(fmax);
xlim([-fmax fmax]*1.2); ylim([-fmax fmax]*1.2);
title('Bartlett-averaged directional spectrum'); colorbar; caxis([0 1]);

figure('Position',[100 100 600 600]);
pcolor(Xc, Yc, E2D_welch_norm); shading interp; colormap(jet); axis equal; hold on;
add_polar_grid(fmax);
xlim([-fmax fmax]*1.2); ylim([-fmax fmax]*1.2);
title('Welch-averaged directional spectrum'); colorbar; caxis([0 1]);

function add_polar_grid(rmax)
    rings = 0:0.05:rmax;
    th = linspace(0,2*pi,200);
    for r = rings
        plot(r*cos(th), r*sin(th), 'Color',[0.7 0.7 0.7],'LineStyle',':');
    end
    for ang = 0:45:315
        plot([0 rmax*sind(ang)], [0 rmax*cosd(ang)], 'Color',[0.7 0.7 0.7],'LineStyle',':');
    end
    text(0,rmax*1.05,'North','HorizontalAlignment','center');
    text(rmax*1.08,0,'East','HorizontalAlignment','left');
    text(0,-rmax*1.1,'South','HorizontalAlignment','center');
    text(-rmax*1.08,0,'West','HorizontalAlignment','right');
end

%% Polar plot

%% Plot E (comparison): Bartlett vs Welch auto-spectra
figure;

plot(freq_axis, Czz_avg, 'LineWidth', 1.2); hold on;
plot(freq_axis, Czz_avg_w, '--', 'LineWidth', 1.2);
xlim([0 0.3]); ylabel('C_{zz} PSD');
legend('Bartlett', 'Welch'); title('Heave (C_{zz})'); grid on;
