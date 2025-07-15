clear;
clc;
close all;

%% --- User parameters ---
filename     = 'id_FRF_20hanning_onlyDisp.txt';   % your file name
headerLines  = 1;                         % skip only the first line

% --- Read all data after the first line ---
% readmatrix will consume everything except the first line
raw = readmatrix(filename, ...
                 'FileType',       'text', ...
                 'Delimiter',      '\t', ...
                 'NumHeaderLines', headerLines);

% keep only rows where the first column is numeric (not NaN)
valid = ~isnan(raw(:,1));
data  = raw(valid, 1:2);

% --- Split into frequency & magnitude arrays ---
freq_data    = data(:,1);     % Hz
mag_lin_data = data(:,2);     % linear magnitude

% --- Convert to dB (ignore non-positive entries) ---
mag_lin_data(mag_lin_data <= 0) = NaN;      % mark zeros/negatives as NaN
mag_dB_data   = 20*log10(mag_lin_data);    % 20·log₁₀ of linear mags

%% --- Step 2: Define your system ---
% Example: a simple first-order low-pass with cutoff at 10 Hz:
%   H(s) = 1/(1 + s/(2*pi*10))
% Replace this with your own numerator/denominator or state-space.
%% Loading Model with Standard Tune

mT=1.9751*1e3; %Platen mass (mp=1.9751 t)
cT=5.78*1e3;   %Total damping, actuator + platen (ct=5.78 kN s/m1)
mass=2e3;
m1 = mass; % kg % 1st mode
f1 = 1.5; % Hz   % 1.5 < f1 < 4
zeta1 = 0.02 ; % 2 < zeta1 < 10
m2 = mass; % kg %2nd mode
f2 =6; % Hz % 6 < f2 < 10
zeta2 = 0.05; % 5 < zeta2 < 25r

c1 = zeta1*2*m1*2*pi*f1; c2 = zeta2*2*m2*2*pi*f2; %N/m/s%coupled 2DOF system

k_p=1.2993/1e-2; %SI units %Pgain (kp=1.2993 V/cm) 
G_c = tf(k_p,1);% Controller

% s,G_T,G_1,G_2,G_T1 ,G_21 ,G_svq,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT , G_Fp_isv  ,c1,c2,k1,k2,AA , BB , CC , DD 
[~,~,~,~,~ ,~ ,~,~,~,~,~,~,sys,~,~ , ~ ,~,~,~,~ ,  ~ ,~,~,~  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);
% Alternatively, if you have transfer function in Hz form:
% sys = tf([gain], [1/(2*pi*f_c) 1]); etc.
% Or define sys via ss(A,B,C,D) as needed.

%% --- Step 3: Choose frequency vector for system response ---
% Cover at least the data range; you may extend a bit beyond for context.
f_min = min(freq_data);
f_max = max(freq_data);
% Create a log-spaced vector of frequencies (in Hz) for a smooth curve:
nPoints = 200;  % increase for smoother curve
f_vec = logspace(log10(f_min*0.8), log10(f_max*1.2), nPoints);  % in Hz
omega_vec = 2*pi * f_vec;  % convert to rad/s

%% --- Step 4: Compute system frequency response ---
% bode(sys, omega) returns magnitude and phase at specified rad/s frequencies.
[mag_sys, phase_sys, wout] = bode(sys, omega_vec);
% mag_sys has size 1×1×nPoints (for single-input single-output). Squeeze to vector:
mag_sys = squeeze(mag_sys);    % linear magnitude
phase_sys = squeeze(phase_sys);% in degrees
% wout is rad/s, should match omega_vec.

% Convert to dB and to Hz:
mag_sys_dB = 20*log10(mag_sys);
freq_sys_Hz = wout / (2*pi);

%% --- Step 5: Plot magnitude overlay ---
figure;
% Plot measured data:
semilogx(freq_data, mag_dB_data, 'o', 'MarkerSize',6, 'LineWidth',1.2);
hold on;
% Plot system curve:
semilogx(freq_sys_Hz, mag_sys_dB, '-', 'LineWidth',1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Measured Data vs. System Bode Magnitude');
legend('Identified in Adapt.exe','Model','Location','best');
xlim([min(f_vec), max(f_vec)]);
% Optionally adjust ylim or let MATLAB autoscale.

%% --- (Optional) Step 6: Plot phase in a second subplot ---
% If you care about phase comparison (but you have no measured phase to overlay),
% you can still visualize system phase:
% figure;
% semilogx(freq_sys_Hz, phase_sys, '-', 'LineWidth',1.5);
% grid on;
% xlabel('Frequency (Hz)');
% ylabel('Phase (degrees)');
% title('System Bode Phase');
% xlim([min(f_vec), max(f_vec)]);

%% --- Alternative: Single figure with 2 subplots (magnitude + phase) ---
%{
figure;
subplot(2,1,1);
semilogx(freq_data, mag_dB_data, 'o', 'MarkerSize',6, 'LineWidth',1.2);
hold on;
semilogx(freq_sys_Hz, mag_sys_dB, '-', 'LineWidth',1.5);
grid on;
ylabel('Magnitude (dB)');
title('Measured vs System Bode');
legend('Measured','System','Location','best');
xlim([min(f_vec), max(f_vec)]);

subplot(2,1,2);
semilogx(freq_sys_Hz, phase_sys, '-', 'LineWidth',1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
title('System Phase');
xlim([min(f_vec), max(f_vec)]);
%}

%% --- Notes & Tips ---
% 1. If your system is multi-input or multi-output, or you want a specific input-output pair,
%    extract the appropriate SISO TF before using bode.
% 2. If you prefer to use MATLAB’s bodeplot with `hold on`, note that bodeplot by default
%    uses rad/s axis. Overlaying measured data on that plot is more cumbersome; 
%    manual semilogx as above is clearer.
% 3. If your measured data has phase as well, you could overlay measured phase points similarly:
%    semilogx(freq_data, phase_data, 'x', ...).
% 4. Ensure frequency units match: `bode` uses rad/s, so convert to Hz as shown.
% 5. If you want to match plotting style, you can customize markers, line styles, colors, etc.
% 6. If your data span is very wide, you may adjust f_vec beyond data range to see roll-off.
% 7. If you want to plot linear-magnitude overlay: plot mag_lin_data vs freq_data with `semilogx`,
%    and plot `mag_sys` vs `freq_sys_Hz` in linear scale (no dB). But Bode convention is dB.
% 8. If your `sys` has time delays or other peculiarities, bode still handles them.
% 9. If you prefer to directly call `bode` and capture the plot handle:
%    h = bodeplot(sys);
%    % then modify axes and overlay data manually by retrieving the axes handle:
%    ax = findall(gcf,'Type','axes');
%    hold(ax(1),'on'); semilogx(ax(1),freq_data, mag_dB_data, 'o');
%    But manual computation as above is simpler.
%
% After you adapt `sys` to your actual transfer function, run the script: you will see your data
% points and the system’s Bode curve overlaid, enabling comparison of measured vs. model.
