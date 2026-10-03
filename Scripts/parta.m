%% pn_junction_part_a.m
% Computes and plots the electric field and electrostatic potential 
% for a single applied voltage (VA = -20 V).

clear; clc; close all;   % for starting fresh every time the script is run

%% 1. Constants & Given Parameters
q    = 1.602e-19;       % Electronic charge (C)
k    = 1.38e-23;        % Boltzmann constant (J/K)
T    = 300;             % Temperature (K)
eps0 = 8.854e-14;       % Vacuum permittivity (F/cm)
epsr = 11.7;            % Relative permittivity of Si
eps  = epsr * eps0;     % Permittivity of Si (F/cm)
ni   = 1.5e10;            % Intrinsic carrier conc. (cm^-3)
Vt   = k * T / q;       % Thermal voltage (V)

NA = 1e15;              % Acceptor doping concentration (cm^-3)
ND = 2e14;              % Donor doping concentration (cm^-3)
VA = -20;               % Applied voltage (V)

%% 2. Depletion Approximation Calculations
Vbi = Vt * log(NA * ND / ni^2);                 % Built-in potential (V)
Vj  = Vbi - VA;                                 % Total potential drop (V)
W   = sqrt(2 * eps * Vj / q * (1/NA + 1/ND));   % Total depletion width (cm)
xn  = W * NA / (NA + ND);                       % n-side width (cm)
xp  = W * ND / (NA + ND);                       % p-side width (cm)
E0  = -q * NA * xp / eps;                       % Electric field at x=0 (V/cm)
V0  = q * NA * xp^2 / (2 * eps);                % Potential at x=0 (V)

%% 3. Position Array and Limits Setup
xmax = 2.5 * max(xn, xp);                       % Maximum x-axis limit (cm)
x = linspace(-xmax, xmax, 2000);                
x = sort([x, 0]);                               % Ensure exactly x=0 is evaluated

%% 4. Profile Generation (Electric Field & Potential)
E = zeros(size(x));
V = zeros(size(x));

% Logical indexing to find points in each region
Lp = (x >= -xp) & (x <= 0);                     % p-side depletion region
Ln = (x > 0) & (x <= xn);                       % n-side depletion region
Rn = (x > xn);                                  % Neutral n-region

% Apply physical formulas to the respective regions
E(Lp) = -q * NA * (x(Lp) + xp) / eps;
E(Ln) = -q * ND * (xn - x(Ln)) / eps;

V(Lp) = q * NA * (x(Lp) + xp).^2 / (2 * eps);
V(Ln) = Vj - q * ND * (xn - x(Ln)).^2 / (2 * eps);
V(Rn) = Vj;                                     % V remains constant past xn

%% 5. Plotting (Subplot coordination)
Emin = 1.1 * E0;                                % y-axis limits per assignment
Vmax = 1.1 * Vj;

figure('Name', 'Part (a) Results', 'Color', 'w');

% Top Plot: Electric Field
subplot(2,1,1);
plot(x * 1e4, E, 'LineWidth', 1.5, 'Color', '#0072BD');
hold on; 
xline(-xp * 1e4, 'r--', 'Alpha', 0.5);
xline(xn * 1e4, 'r--', 'Alpha', 0.5);
plot(-xp * 1e4, 0, 'ko', 'MarkerFaceColor', 'k'); 
text(-xp * 1e4, 0, sprintf('  (%.2f, 0)', -xp * 1e4), 'VerticalAlignment', 'bottom');
plot(xn * 1e4, 0, 'ko', 'MarkerFaceColor', 'k');
text(xn * 1e4, 0, sprintf('  (%.2f, 0)', xn * 1e4), 'VerticalAlignment', 'bottom');
plot(0, E0, 'ko', 'MarkerFaceColor', 'k');
text(0, E0, sprintf('  (0, %.0f)', E0), 'VerticalAlignment', 'top');

xlim([-xmax xmax] * 1e4);
ylim([Emin 0]);
grid on; box on;
xlabel('Position x (\mum)');
ylabel('Electric field E (V/cm)');
title(sprintf('Electric field v/s position: N_A = 10^{15}/cm^3, N_D = 2\\times10^{14}/cm^3, V_A = %g V', VA));

% Bottom Plot: Electrostatic Potential
subplot(2,1,2);
plot(x * 1e4, V, 'LineWidth', 1.5, 'Color', '#D95319');
hold on; 
xline(-xp * 1e4, 'r--', 'Alpha', 0.5);
xline(xn * 1e4, 'r--', 'Alpha', 0.5);
plot(-xp * 1e4, 0, 'ko', 'MarkerFaceColor', 'k');
text(-xp * 1e4, 0, sprintf('  (%.2f, 0)', -xp * 1e4), 'VerticalAlignment', 'bottom');
plot(0, V0, 'ko', 'MarkerFaceColor', 'k');
text(0, V0, sprintf(' (0, %.2f)  ', V0), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
plot(xn * 1e4, Vj, 'ko', 'MarkerFaceColor', 'k');
text(xn * 1e4, Vj, sprintf('  (%.2f, %.2f)', xn * 1e4, Vj), 'VerticalAlignment', 'top');

xlim([-xmax xmax] * 1e4);
ylim([0 Vmax]);
grid on; box on;
xlabel('Position x (\mum)');
ylabel('Potential V (V)');
title(sprintf('Potential v/s Position: N_A = 10^{15}/cm^3, N_D = 2\\times10^{14}/cm^3, V_A = %g V', VA));