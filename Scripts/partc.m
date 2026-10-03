%% partc.m
% Computes and plots the electric field and electrostatic potential 
% for multiple applied voltages and outputs a tabulated list of parameters.

clear; clc; close all;

%% 1. Constants & Base Parameters
q    = 1.602e-19;       
k    = 1.38e-23;        
T    = 300;             
eps0 = 8.854e-14;       
epsr = 11.7;            
eps  = epsr * eps0;     
ni   = 1e10;            
Vt   = k * T / q;       

NA  = 1e15;              
ND  = 2e14;              
VA0 = -20;               

%% 2. Setup Figure and Baseline Limits
Vbi = Vt * log(NA * ND / ni^2);

% Calculate max limits based on n=0 (the -20V case creates the widest boundaries)
Vj_max = Vbi - VA0;
W_max  = sqrt(2 * eps * Vj_max / q * (1/NA + 1/ND));
xn_max = W_max * NA / (NA + ND);
xp_max = W_max * ND / (NA + ND);
xmax   = 2.5 * max(xn_max, xp_max);
E0_max = -q * NA * xp_max / eps;

x = linspace(-xmax, xmax, 2000);                
x = sort([x, 0]);

figure('Name', 'Part (c) Results', 'Color', 'w');
legend_labels = cell(1, 4);

%% 3. Print Header & Loop Through Voltages
% Dynamically format NA and ND to look like '1x10^15' instead of '1.00e+15'
str_NA = sprintf('%gx10^%d', NA/10^floor(log10(NA)), floor(log10(NA)));
str_ND = sprintf('%gx10^%d', ND/10^floor(log10(ND)), floor(log10(ND)));

% Create the table header in the Command Window
fprintf('\n================================================ PART (c) OUTPUT ================================================\n');
fprintf('  V_A (V) | N_A (cm^-3) | N_D (cm^-3) | V_bi (V) | x_n (um) | x_p (um) |  W (um)  | E at x=0 (V/cm)| V at x=0 (V)\n');
fprintf('-----------------------------------------------------------------------------------------------------------------\n');

for n = 0:3
    VA = VA0 / (2^n);               
    legend_labels{n+1} = sprintf('V_A = %.2f V', VA);

    % Recalculate parameters for this specific VA
    Vj = Vbi - VA;
    W  = sqrt(2 * eps * Vj / q * (1/NA + 1/ND));
    xn = W * NA / (NA + ND);
    xp = W * ND / (NA + ND);

    E0 = -q * NA * xp / eps;
    V0 = q * NA * xp^2 / (2 * eps);

    % Print the requested parameters for this iteration to the Command Window
    % Multiply xn, xp, and W by 1e4 to display them clearly in micrometers (um)
    fprintf('%9.2f | %11s | %11s | %8.4f | %8.3f | %8.3f | %8.3f | %14.0f | %12.4f\n', ...
        VA, str_NA, str_ND, Vbi, xn*1e4, xp*1e4, W*1e4, E0, V0);

    E = zeros(size(x));
    V = zeros(size(x));

    Lp = (x >= -xp) & (x <= 0);
    Ln = (x > 0) & (x <= xn);
    Rn = (x > xn);

    E(Lp) = -q * NA * (x(Lp) + xp) / eps;
    E(Ln) = -q * ND * (xn - x(Ln)) / eps;

    V(Lp) = q * NA * (x(Lp) + xp).^2 / (2 * eps);
    V(Ln) = Vj - q * ND * (xn - x(Ln)).^2 / (2 * eps);
    V(Rn) = Vj;

    % --- Plot Electric Field ---
    subplot(2,1,1);
    hold on; 
    pE = plot(x * 1e4, E, 'LineWidth', 1.5);
    c = pE.Color; 

    plot([-xp*1e4, xn*1e4], [0, 0], 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(0, E0, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerSize', 5, 'HandleVisibility', 'off');

    % --- Plot Potential ---
    subplot(2,1,2);
    hold on; 
    plot(x * 1e4, V, 'LineWidth', 1.5, 'Color', c);

    plot(-xp*1e4, 0, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(0, V0, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(xn*1e4, Vj, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerSize', 5, 'HandleVisibility', 'off');
end

%% 4. Apply Final Formatting
subplot(2,1,1);
xlim([-xmax xmax] * 1e4);
ylim([1.1*E0_max 0]);
grid on; box on;
xlabel('Position x (\mum)');
ylabel('Electric field E (V/cm)');
title(sprintf('Si pn step junction, T = 300 K: N_A = 10^{15} cm^{-3}, N_D = 2\\times10^{14} cm^{-3}'));
legend(legend_labels, 'Location', 'best');

subplot(2,1,2);
xlim([-xmax xmax] * 1e4);
ylim([0 1.1*Vj_max]);
grid on; box on;
xlabel('Position x (\mum)');
ylabel('Potential V (V)');
legend(legend_labels, 'Location', 'best');