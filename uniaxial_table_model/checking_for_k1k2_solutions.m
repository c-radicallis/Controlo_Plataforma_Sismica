clear; clc; close all;

% Define constants (masses in kg)
m1 = 2000; 
m2 = 2000; 

% Define lists of frequencies (Hz)
f1_list = [1.5, 2 , 3 , 4 ];  
f2_list = 2.415*f1_list;%[ 6,  8 , 10];  % Ensure these vectors have the same length

% Define range for k1 and k2 (stiffness in N/m)
k1_vals = linspace(1e3, 8e6, 100);
k2_vals = linspace(1e3, 4e6, 100);
[K1, K2] = meshgrid(k1_vals, k2_vals);

figure(1);
hold on;
% Use a colormap to assign different colors
colors = lines(length(f1_list)*2); % Two colors per pair (one for each contour)

for idx = 1:length(f1_list)
    f1 = f1_list(idx);
    f2 = f2_list(idx);
    
    % Compute squared circular frequencies for each f
    omega1_sq = (2*pi*f1)^2;
    omega2_sq = (2*pi*f2)^2;
    
    % Common term for both equations
    term1 = (m1 .* K2 + m2 .* (K1 + K2)) / (2*m1*m2);
    term2 = 0.5 * sqrt(((m1 .* K2 + m2 .* (K1 + K2))/(m1*m2)).^2 - 4*(K1 .* K2)/(m1*m2));
    
    % Define equations (set equal to zero)
    eq1 = term1 - term2 - omega1_sq;
    eq2 = term1 + term2 - omega2_sq;
    
    % Plot contour for eq1 (corresponding to f1)
    contour(K1, K2, eq1, [0 0], 'LineColor', colors(2*idx-1,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('f1 = %.1f Hz', f1));
    
    % Plot contour for eq2 (corresponding to f2)
    contour(K1, K2, eq2, [0 0],"--", 'LineColor', colors(2*idx,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('f2 = %.1f Hz', f2));
end

xlabel('Stiffness k_1 (N/m)');
ylabel('Stiffness k_2 (N/m)');
title('Contours for Multiple f_1 and f_2 Values');
legend('show');
grid on;
hold off;