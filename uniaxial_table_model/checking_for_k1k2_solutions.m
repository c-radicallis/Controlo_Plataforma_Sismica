clear; clc; close all;

%% Define constants
m1 = 2000; % kg
m2 = 2000; % kg

%% Define ranges for k1 and k2 and create a grid
k1_vals = linspace(1e3, 1e7, 1000);  % Adjust range and resolution as needed
k2_vals = linspace(1e3, 1e7, 1000);
[K1, K2] = meshgrid(k1_vals, k2_vals);

%% Pre-compute terms common to both equations
term1 = (m1.*K2 + m2.*(K1+K2)) / (2*m1*m2);
% Compute the discriminant inside the square root:
sqterm = ((m1.*K2 + m2.*(K1+K2))/(m1*m2)).^2 - 4*(K1.*K2)/(m1*m2);
% Avoid complex results by setting negative values to NaN:
sqterm(sqterm < 0) = NaN;
term2 = 0.5 * sqrt(sqterm);

%% Define list of [f1, f2] pairs (Hz)
% Each row: [f1, f2]
freqPairs = [2.5, 6; 3,7 ; 4,8 ; 5,10];  % Add more pairs as desired

%% Create figure and loop over frequency pairs to plot contours
figure;
hold on;
colors = lines(size(freqPairs,1));  % Get distinct colors for each pair

for idx = 1:size(freqPairs,1)
    % Extract frequency pair for this iteration
    f1 = freqPairs(idx,1);
    f2 = freqPairs(idx,2);
    
    % Compute squared angular frequencies
    omega1_sq = (2*pi*f1)^2;
    omega2_sq = (2*pi*f2)^2;
    
    % Define the two equations:
    % Equation for f1: term1 - term2 = (2*pi*f1)^2
    % Equation for f2: term1 + term2 = (2*pi*f2)^2
    eq1 = term1 - term2 - omega1_sq;
    eq2 = term1 + term2 - omega2_sq;
    
    % Plot the zero-level contour for eq1 (solid line)
    [~, h1] = contour(K1, K2, eq1, [0 0], 'LineColor', colors(idx,:), 'LineWidth', 2);
    h1.DisplayName = sprintf('f1 = %.2f Hz', f1);
    
    % Plot the zero-level contour for eq2 (dashed line)
    [~, h2] = contour(K1, K2, eq2, [0 0], '--', 'LineColor', colors(idx,:), 'LineWidth', 2);
    h2.DisplayName = sprintf('f2 = %.2f Hz', f2);
end

xlabel('Stiffness k_1 (N/m)');
ylabel('Stiffness k_2 (N/m)');
title('Contour Lines for Multiple f_1 and f_2 Pairs');
legend('show');
grid on;
hold off;
