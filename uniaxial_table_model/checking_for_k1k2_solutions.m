clear; clc; close all;

% Define given constants
m1 = 2000; % kg
m2 = 2000; % kg
f1 = 2.4;  % Hz
f2 = 6;    % Hz

% Define frequency terms squared
omega1_sq = (2 * pi * f1)^2;
omega2_sq = (2 * pi * f2)^2;

% Define range for k1 and k2
k1_vals = linspace(1e3, 1e7, 1000);  % Stiffness range for k1
k2_vals = linspace(1e3, 1e7, 1000);  % Stiffness range for k2

% Create meshgrid
[K1, K2] = meshgrid(k1_vals, k2_vals);

% Compute left-hand side of equations
term1 = (m1 .* K2 + m2 .* (K1 + K2)) ./ (2 * m1 * m2);
term2 = sqrt(((m1 .* K2 + m2 .* (K1 + K2)) ./ (m1 * m2)).^2 - 4 * (K1 .* K2) ./ (m1 * m2)) / 2;

% Define the two equations
eq1 = term1 - term2 - omega1_sq;
eq2 = term1 + term2 - omega2_sq;

% Plot contour where each equation is zero
figure(1);
hold on;
contour(K1, K2, eq1, [0 0], 'LineWidth', 1, 'DisplayName', sprintf('Equation 1 - %.1f %.1f' , f1, f2));
contour(K1, K2, eq2, [0 0], 'LineWidth', 1, 'DisplayName', sprintf('Equation 2 - %.1f %.1f' , f1, f2));
xlabel('Stiffness k_1 (N/m)');
ylabel('Stiffness k_2 (N/m)');
title('Solution Curves for k_1 and k_2');
legend;
grid on;
hold off;
