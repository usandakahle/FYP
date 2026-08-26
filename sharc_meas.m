%% SHARC buoy measurement 
% Takes the  open-ocean surface elevation eta(t) generated in
% ocean_wave_spectrum.m and passes it through the SHARC buoy's actual
% sensing + onboard DSP chain, so you can see what the buoy would have
% reported (Hm0, Tm01, PSD) versus the ground-truth input spectrum.
%
% Only heave (z-axis accel) is used below, matching the current SHARC
% pipeline — surge/sway are captured by the IMU but not processed into
% the wave-height product (that's the direction problem I still have
% to add).
%
%
% Requires the Signal Processing Toolbox (fir1, filtfilt, hanning).

run('ocean_wave_spectrum.m');   % brings eta, t, f, a, alpha, fs, Hs, Tp into scope


%% Heave Acceleration
% eta(t) = sum_i a(i)*cos(2*pi*f(i)*t + alpha(i))
% Differentiate the same sum analytically (avoids numerical-diff noise):
% d^2/dt^2 [a*cos(wt+phi)] = -a*w^2*cos(wt+phi)
w = 2*pi*f;
fs_imu = 100;                     % (Hz) ICM20649 sample rate
t_imu = 0: 1/fs_imu: D-1/fs_imu;

heave_spec = zeros(size(t_imu));
for i = 1:length(f)
    heave_spec = heave_spec - a(i)*w(i)^2*cos(w(i)*t_imu + alpha(i));
end

%% IMU sensor model (ICM20649)


noise_density = 289e-6 * 9.81;   % 289 ug/sqrt(Hz) -> m/s^2/sqrt(Hz)
sigma_noise = noise_density * sqrt(fs_imu/2);
accel_noise = sigma_noise * randn(size(t_imu));
% 32-bit stored, 4g full scale: quantisation step is negligible at these
% amplitudes (waves-in-ice < 1 m/s^2), so it's omitted here.

heave = heave_spec + accel_noise;

%% Already Existing DSP pipeline 
% Detrend raw block(removes mean data)
heave= detrend(heave);

% FIR low-pass + decimate, M = 32 (100 Hz -> 3.125 Hz nominal)
M = 32;
fs_dec = fs_imu / M;
if fs_imu/2 <= fs_dec
    warning(['fs_imu = %.2f Hz cannot be meaningfully decimated by M=%d ' ...
        '(Nyquist is already at/under the target rate). Numbers below ' ...
        'are structural only — see the mismatch note at the top.'], fs_imu, M);
end
cutoff_norm = min(0.9, 1.12 / (fs_imu/2));   % thesis cutoff = 1.12 Hz
lpf = fir1(64, cutoff_norm);
heave_filt = filtfilt(lpf, 1, heave);
heave_dec = downsample(heave_filt, M);

%  RC high-pass detrend, tau = 100 s (thesis: k=0.9986 @ 6.25 Hz;
% same formula generalised to whatever fs_dec is found t \o be
tau = 100;
k = exp(-1/(tau*fs_dec));
heave_hp = zeros(size(heave_dec));
for n = 2:length(heave_dec)
    heave_hp(n) = k*heave_hp(n-1) + k*(heave_dec(n) - heave_dec(n-1));
end

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

win = hanning(Nfft);                   %windowing function
%win_energy_correction = 1.63;      % thesis-quoted Hanning energy correction
freq_axis = (0:Nfft/2-1) * (fs_dec/Nfft);
%freq_axis(1) = freq_axis(2)/2;     % avoid divide-by-zero at DC below

seg_psds = zeros(Nfft/2, Nseg);
for s = 1:Nseg
    idx = (s-1)*Nfft + (1:Nfft);
    x = heave_hp(idx) .* win(:)';

    X = fft(x, Nfft);
    X = X(1:Nfft/2);                       % one-sided

    % acceleration -> displacement: divide by -(2*pi*f)^2
    valid = freq_axis >= f_min & freq_axis <= 0.5;
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
    f_min = 0.03;   % Hz, thesis drift cutoff (~33 s period)
    X_disp = X_disp .* (freq_axis >= f_min & freq_axis <= 0.5);

    %  PSD, energy-corrected, Parseval-scaled
    Pxx = (2 / (fs_dec * sum(win.^2))) * abs(X_disp).^2;
    
    seg_psds(:,s) = Pxx;
end

%%  Bartlett-averaged PSD (consider welch)
Pxx_avg = mean(seg_psds, 2);


%%  Wave parameters from the SHARC-measured PSD 
m0 = trapz(freq_axis, Pxx_avg);
m1 = trapz(freq_axis, freq_axis .* Pxx_avg);
Hm0_sharc = 4*sqrt(m0);
Tm01_sharc = m0./m1;

fprintf('SHARC-measured : Hm0 = %.3f m, Tm01 = %.2f s\n', Hm0_sharc, Tm01_sharc);
fprintf('Input truth    : Hs  = %.3f m, Tp  = %.2f s (peak, not Tm01 — compare with care)\n', Hs, Tp);

%%  Plots
figure;

%plot(t, heave_spec, t, heave, '--');
%xlabel('Time (s)'); ylabel('Heave accel (m/s^2)');
%legend('true heave','IMU-measured heave');
%title('Heave acceleration: truth vs sensor');

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

figure;
plot(freq_axis, Pxx_avg);
xlabel('Frequency (Hz)'); 
ylabel('Displacement PSD (m^2/Hz)');
title('SHARC-recovered heave displacement spectrum');