%% direction validation plots

run('ocean_wave_spectrum.m');   % brings f, a, alpha, theta_i, D_theta, theta_grid,
                                 % theta_mean, sigma, s, eta, Sx, Sy, t, fs, E_f into scope

%% Plot 1: Target directional spectrum E(f,theta) — compass-bearing heatmap
theta_compass_deg      = mod(90 - rad2deg(theta_grid), 360);
theta_mean_compass_deg = mod(90 - rad2deg(theta_mean), 360);

E_f_theta = E_f(:) .* D_theta(:)';   % outer product: E(f) x D(theta)

fmax   = 0.3;
band   = f <= fmax;
f_band = f(band);
E_norm = E_f_theta(band,:) / max(E_f_theta(band,:), [], 'all');

[Theta, F] = meshgrid(theta_compass_deg, f_band);
Xc = F .* sind(Theta);
Yc = F .* cosd(Theta);

figure('Position',[100 100 600 600]);
pcolor(Xc, Yc, E_norm); shading interp; colormap(jet); axis equal; hold on;
plot([0 fmax*sind(theta_mean_compass_deg)], [0 fmax*cosd(theta_mean_compass_deg)], ...
    'k--', 'LineWidth', 1.5);
add_polar_grid(fmax);
xlim([-fmax fmax]*1.2); ylim([-fmax fmax]*1.2);
title(sprintf('Target directional spectrum E(f,\\theta), s = %.2f', s));
colorbar; caxis([0 1]);

%% Plot 2: Heave vs slope quadrature relationship (short window)
idx_win = t <= 100;
figure;
plot(t(idx_win), eta(idx_win), 'LineWidth', 1.2); hold on;
plot(t(idx_win), Sx(idx_win), 'LineWidth', 1.0);
plot(t(idx_win), Sy(idx_win), 'LineWidth', 1.0);
xlabel('Time (s)');
ylabel('Heave (m) / Slope (rad)');
legend('\eta(t) heave', 'S_x(t) slope', 'S_y(t) slope');
title('Heave and surface slope components (first 100 s)');
grid on;

%% Plot 3: Energy-weighted directional histogram vs target spreading (compass heatmap-style polar)
nbins = 36;
edges = linspace(-pi, pi, nbins+1);
centers = edges(1:end-1) + diff(edges)/2;

bin_idx = discretize(theta_i, edges);
energy_per_bin = accumarray(bin_idx(:), (a.^2)', [nbins, 1]);
energy_per_bin = energy_per_bin / trapz(centers, energy_per_bin);

centers_compass_deg = mod(90 - rad2deg(centers), 360);
[~, sort_idx] = sort(centers_compass_deg);
centers_compass_sorted = centers_compass_deg(sort_idx);
energy_sorted = energy_per_bin(sort_idx);

figure;
pax = polaraxes;
polarplot(pax, deg2rad([centers_compass_sorted centers_compass_sorted(1)]), ...
    [energy_sorted; energy_sorted(1)], 'LineWidth', 1.5);
hold on;
polarplot(pax, deg2rad(theta_compass_deg), D_theta, 'r--', 'LineWidth', 1.5);
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
title('Energy-weighted draws of \theta_i vs target D(\theta) (compass bearing)');
legend('Drawn directions (energy-weighted)', 'Target D(\theta)', 'Location', 'southoutside');

%% Plot 4: Heave PSD vs target E(f) (sanity check)
[Pxx_eta, f_est] = pwelch(eta, hamming(1024), 512, [], fs);
figure;
plot(f, E_f, 'LineWidth', 1.5); hold on;
plot(f_est, Pxx_eta, 'LineWidth', 1.2);
xlim([0 1]);
xlabel('Frequency (Hz)'); ylabel('S(f) (m^2/Hz)');
legend('Target E(f)', 'Recovered heave PSD');
title('Heave spectrum check');
grid on;

%% Plot 5: Slope PSDs vs theoretical prediction
cos2_avg = trapz(theta_grid, D_theta .* cos(theta_grid - theta_mean).^2);
sin2_avg = trapz(theta_grid, D_theta .* sin(theta_grid - theta_mean).^2);
k_of_f = (2*pi*f).^2 / g;
E_Sx_theory = E_f .* k_of_f.^2 * cos2_avg;
E_Sy_theory = E_f .* k_of_f.^2 * sin2_avg;

[Pxx_Sx, ~] = pwelch(Sx, hamming(1024), 512, [], fs);
[Pxx_Sy, ~] = pwelch(Sy, hamming(1024), 512, [], fs);

figure;
subplot(2,1,1);
plot(f, E_Sx_theory, 'LineWidth', 1.5); hold on;
plot(f_est, Pxx_Sx, 'LineWidth', 1.2);
xlim([0 1]); xlabel('Frequency (Hz)'); ylabel('S_x PSD (rad^2/Hz)');
legend('Theoretical', 'Recovered'); title('Slope S_x spectrum check'); grid on;

subplot(2,1,2);
plot(f, E_Sy_theory, 'LineWidth', 1.5); hold on;
plot(f_est, Pxx_Sy, 'LineWidth', 1.2);
xlim([0 1]); xlabel('Frequency (Hz)'); ylabel('S_y PSD (rad^2/Hz)');
legend('Theoretical', 'Recovered'); title('Slope S_y spectrum check'); grid on;

%% Plot 6: Direction recovery self-check (heave-slope cross-spectra)
nfft = 1024; noverlap = 512; win = hamming(nfft);
[Czz, fcp] = cpsd(eta, eta, win, noverlap, [], fs);
[Cxx, ~]   = cpsd(Sx, Sx, win, noverlap, [], fs);
[Cyy, ~]   = cpsd(Sy, Sy, win, noverlap, [], fs);
[Czx, ~]   = cpsd(eta, Sx, win, noverlap, [], fs);
[Czy, ~]   = cpsd(eta, Sy, win, noverlap, [], fs);
[Cxy, ~]   = cpsd(Sx, Sy, win, noverlap, [], fs);

Qzx = imag(Czx);
Qzy = imag(Czy);

a1 = Qzx ./ sqrt(Czz .* (Cxx + Cyy));
b1 = Qzy ./ sqrt(Czz .* (Cxx + Cyy));
theta_recovered = atan2(b1, a1);

figure;
plot(fcp, rad2deg(wrapToPi(theta_recovered)), 'LineWidth', 1.2); hold on;
yline(rad2deg(theta_mean), 'r--', 'LineWidth', 1.5);
xlim([0.02 0.3]);
xlabel('Frequency (Hz)'); ylabel('Recovered mean direction (deg)');
legend('\theta_{recovered}(f)', '\theta_{mean} (input)');
title('Direction recovery self-check (clean synthetic signal, no IMU/DSP noise)');
grid on;

%% ---- Local functions (must be at end of script) ----
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