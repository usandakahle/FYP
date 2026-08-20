%Generating simulation spectrum from Holthujsen description
%of oceanic waves

lat = ncread('Era5.nc','latitude');
lon = ncread('Era5.nc','longitude');
time = ncread('Era5.nc','valid_time');
swh  = ncread('Era5.nc','swh');   % significant wave height
pp1d = ncread('Era5.nc','pp1d');  % peak wave period

%Locating coordinates in the era5 file (using Ardhuin coordinates)
target_lat  = -60;                          % degrees N (negative = South)
target_lon  = 355;                          % degrees E, 0-360 convention (355 = 5W)
target_time = datenum(2023,6,15) ...        % example date, austral winter
              *86400 - datenum(1970,1,1)*86400;  % convert to seconds since 1970-01-01, matching valid_time units


%Simulation parameters
D = 1000;         %Total duration of record (seconds)
fs = 2;           %sampling frequency (Hz)
dt = 1/fs;        %Time step (seconds)
t = 0:dt:D-dt;    %Time vector
N_t = length(t);

%Frequency Discretization
df = 1/D;         %Frequency interval (Hz) eq.3.5.1
f_max = fs/2;     % Nyquist frequency
f = df:df:f_max;  %frequency vector
N_f = length(f);

%Target variance density spectrum E(f)
%Pierson-Moskowitz Spectrum (Change it for the 
%a spectrum matching the antarctic miz or similar)

% find nearest grid point/time to buoy location & event


[~, ilat] = min(abs(lat - target_lat));
[~, ilon] = min(abs(lon - target_lon));
swh_series = squeeze(swh(ilon, ilat, :));   % Hs over your whole downloaded window

[max_hs, itime] = max(swh_series);          % find the peak storm within the window
target_time = time(itime);                  % valid_time is already in the right units



Hs = swh(ilon, ilat, itime);      %Significant wave height (m)
Tp = pp1d(ilon, ilat, itime);       %Peak period(s)
fp = 1/Tp;        %Peak frequency (Hz)


%E_pm(f) formula
E_f = (0.3125 * Hs^2* fp^4 ./ f.^5).*exp(-1.25 * (fp./f).^4);

% Random-phase/ampliutude model
alpha = 2*pi* rand(1,N_f);

%Rayleigh distribution
b = sqrt(E_f *df);
a = raylrnd(repmat(b,1,1));     %generates one realtization of amplitudes

%Summation to obtain surface elevation
eta = zeros(size(t));
for i = 1:N_f
    eta = eta + a(i) * cos(2*pi*f(i)*t + alpha(i));
end

%plots
figure;
subplot(2,1,1);
plot(f,E_f,'LineWidth',1.5);
xlabel('Frequency (Hz)');
ylabel('Spectral Density E(f) (m^2/Hz)');
title('Target Variance density Spectrum');
grid on;

subplot(2,1,2);
plot(t, eta, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Surface Elevation (m)');
title('Simulated Ocean Surface Elevation');
grid on;

fprintf('Signifcant wave height: %.2f m\n', 4*std(eta));

%Validate eta against E_f using welch (Windowing) reducing noise
[Pxx, f_est] = pwelch(eta, hamming(1024), 512, [], fs);
%

figure;
plot(f, E_f, 'LineWidth', 1.5); hold on;
plot(f_est, Pxx, 'LineWidth', 1.2);
legend('Target E(f)', 'Recovered spectrum (Welch)');
xlabel('Frequency (Hz)'); ylabel('Spectral Density (m^2/Hz)');
title('Target vs. Recovered Spectrum'); grid on;


%Important parameters for the directional calculations.
%Spectral moment sdervied from periods 
m0 = trapz(f, E_f);
m1 = trapz(f, f.*E_f);
m2 = trapz(f, f.^2 .*E_f);
m4 = trapz(f, f.^4 .*E_f);

Hs_spec = 4*sqrt(m0);
Tm01 = m0/m1;
Tm02 = sqrt(m0/m2);

%spectral bandwidth narrowness parameter
nu = sqrt(m0*m2/m1^2 - 1);   % Longuet-Higgins spectral width parameter

%Zero crossing analysis for individual wave heights
zero_cross = find(eta(1:end-1) < 0 & eta(2:end) >= 0);
wave_heights = zeros(length(zero_cross)-1,1);
for k = 1:length(zero_cross)-1
    seg = eta(zero_cross(k):zero_cross(k+1));
    wave_heights(k) = max(seg) - min(seg);
end

figure;
histogram(wave_heights, 'Normalization','pdf'); hold on;
Hrms = sqrt(mean(wave_heights.^2));
H_range = linspace(0, max(wave_heights), 100);
plot(H_range, (2*H_range/Hrms^2).*exp(-(H_range/Hrms).^2));
legend('Empirical', 'Theoretical Rayleigh');
title('Individual Wave Height Distribution');