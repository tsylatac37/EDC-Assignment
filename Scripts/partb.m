%% pn_junction_part_b.m
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

figure('Name', 'Part (b) Results', 'Color', 'w');

% Arrays to store legend labels
legend_labels = cell(1, 4);

%% 3. Loop Through Voltages (n = 0, 1, 2, 3)
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
    c = pE.Color; % Extract the auto-assigned color for this specific line

    % Draw dots but hide them from the legend
    plot([-xp*1e4, xn*1e4], [0, 0], 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(0, E0, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerSize', 5, 'HandleVisibility', 'off');

    % --- Plot Potential ---
    subplot(2,1,2);
    hold on; 
    % Force the potential line to use the exact same color 'c' as the top plot
    plot(x * 1e4, V, 'LineWidth', 1.5, 'Color', c);

    % Draw dots but hide them from the legend
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
title(sprintf('Si pn step junction, Electric field v/s Position T = 300 K: N_A = 10^{15} cm^{-3}, N_D = 2\\times10^{14} cm^{-3}'));
legend(legend_labels, 'Location', 'best');

subplot(2,1,2);
xlim([-xmax xmax] * 1e4);
ylim([0 1.1*Vj_max]);
grid on; box on;
xlabel('Position x (\mum)');
ylabel('Potential V (V)');
title(sprintf('Si pn step junction, Potential v/s Position T = 300 K: N_A = 10^{15} cm^{-3}, N_D = 2\\times10^{14} cm^{-3}'));
legend(legend_labels, 'Location', 'best');