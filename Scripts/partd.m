%% partd.m
% Examines how electric field and potential vary with relative doping magnitudes.

clear; clc; close all;

%% 1. Constants
q    = 1.602e-19;       
k    = 1.38e-23;        
T    = 300;             
eps0 = 8.854e-14;       
epsr = 11.7;            
eps  = epsr * eps0;     
ni   = 1e10;            
Vt   = k * T / q;       
VA   = -20; % Keeping VA constant to isolate doping effects

%% 2. Doping Combinations (NA, ND)
% From assignment prompt part (d)
doping_pairs = [
    1e17, 1e15;
    1e16, 1e15;
    1e15, 1e15;
    1e15, 1e16;
    1e15, 1e17;
    1e18, 1e16;
    1e16, 1e14
];
num_pairs = size(doping_pairs, 1);

% Pre-calculate maximum x-limits to ensure all plots fit on one scale
xmax = 0;
E0_min = 0; 
for i = 1:num_pairs
    NA = doping_pairs(i,1);
    ND = doping_pairs(i,2);
    Vbi = Vt * log(NA * ND / ni^2);
    Vj  = Vbi - VA;
    W   = sqrt(2 * eps * Vj / q * (1/NA + 1/ND));
    xn  = W * NA / (NA + ND);
    xp  = W * ND / (NA + ND);
    E0  = -q * NA * xp / eps;
    
    xmax = max(xmax, 2.5 * max(xn, xp));
    E0_min = min(E0_min, E0);
end

% Use high resolution array to capture extremely narrow depletion regions
x = linspace(-xmax, xmax, 10000); 
x = sort([x, 0]);

figure('Name', 'Part (d) Results', 'Color', 'w');
legend_labels = cell(1, num_pairs);
colors = lines(num_pairs);

%% 3. Calculate and Plot Profiles
for i = 1:num_pairs
    NA = doping_pairs(i,1);
    ND = doping_pairs(i,2);
    
    % Format legend text nicely (e.g., 1e15 instead of 1000000000000000)
    legend_labels{i} = sprintf('N_A = 10^{%d}, N_D = 10^{%d}', log10(NA), log10(ND));
    
    Vbi = Vt * log(NA * ND / ni^2);
    Vj  = Vbi - VA;
    W   = sqrt(2 * eps * Vj / q * (1/NA + 1/ND));
    xn  = W * NA / (NA + ND);
    xp  = W * ND / (NA + ND);
    
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
    
    subplot(2,1,1); hold on;
    plot(x * 1e4, E, 'LineWidth', 1.5, 'Color', colors(i,:));
    
    subplot(2,1,2); hold on;
    plot(x * 1e4, V, 'LineWidth', 1.5, 'Color', colors(i,:));
end

%% 4. Formatting
subplot(2,1,1);
xlim([-xmax xmax] * 1e4);
ylim([1.1*E0_min 0]);
grid on; box on;
xlabel('Position x (\mum)');
ylabel('Electric field E (V/cm)');
title('Effect of Asymmetric Doping on Electric Field (V_A = -20 V)');
% Place legend outside so it doesn't cover the dense curves
legend(legend_labels, 'Location', 'eastoutside'); 

subplot(2,1,2);
xlim([-xmax xmax] * 1e4);
ylim([0 22]); % Max potential drop across all cases is slightly over 21V
grid on; box on;
xlabel('Position x (\mum)');
ylabel('Potential V (V)');
title('Effect of Asymmetric Doping on Potential');
legend(legend_labels, 'Location', 'eastoutside');