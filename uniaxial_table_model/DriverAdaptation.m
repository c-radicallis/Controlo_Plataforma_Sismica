clear;
clc;
close all;

%%
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'

mT=1.9751*1e3; %Platen mass (mp=1.9751 t)
cT=5.78*1e3;   %Total damping, actuator + platen (ct=5.78 kN s/m1)
mass=2e3;
m1 = mass; % kg % 1st mode
f1 = 4; % Hz   % 1.5 < f1 < 4
zeta1 = 0.02 ; % 2 < zeta1 < 10
m2 = mass; % kg %2nd mode
f2 =10; % Hz % 6 < f2 < 10
zeta2 = 0.05; % 5 < zeta2 < 25r

k_p=1.2993/1e-2; %SI units %Pgain (kp=1.2993 V/cm) 
G_c = tf(k_p,1);% Controller

% s,G_T,G_1,G_2,G_T1 ,G_21 ,G_svq,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT , G_Fp_isv  ,c1,c2,k1,k2, AA,BB,CC,DD
[s,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ , ~,~,~,~  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);

%%  Load seismic signal and scale down if necessary
dados = load('elcentro.txt');
t_vector = dados(:,1);
t_step = t_vector(2);
ddx_ref = dados(:,2);

lim_displacement = 0.1; % m % Limits
lim_velocity = 0.4; % m/s
lim_force = 200e3; % N

s=tf('s') ;
x_ref = lsim(1/s^2,  ddx_ref , t_vector ,'foh');
max_xref = max(x_ref);
scale=1;
while max_xref > lim_displacement % Scaling down if necessary
    scale = 0.95*scale;
    ddx_ref = 0.95*ddx_ref;
    x_ref = lsim(1/s^2,  ddx_ref , t_vector ,'foh');
    max_xref = max(x_ref);
end
v_ref =  lsim(1/s,  ddx_ref , t_vector ,'foh');
max_vref = max(v_ref);

%% Plots
close all;
fig2 = figure(2); ax2 = axes(fig2); hold(ax2, 'on'); grid on; hold on; title('Input to Servo');  legend(); xlabel('Time (s)');  ylabel('Voltage (V)');
fig3 = figure(3); ax3 = axes(fig3); hold(ax3, 'on'); grid on; title(' Platen Displacement ');  hold on; legend(); xlabel('Time (s)');  ylabel('Displacement (mm)');
fig4 = figure(4); ax4 = axes(fig4); hold(ax4, 'on'); grid on; title('Platen Acceleration');  hold on; legend(); xlabel('Time (s)');  ylabel('Acceleration (m/s^2)'); 
fig5 = figure(5); ax5 = axes(fig5); hold(ax5, 'on');grid on;title(' Platen Displacement Tracking Error '); legend();xlabel('Time (s)'); ylabel('Error (mm)');
fig6 = figure(6); ax6 = axes(fig6); hold(ax6, 'on');grid on;title(' Platen Acceleration Tracking Error '); legend();xlabel('Time (s)'); ylabel('Error (m/s^2)');
fig7 = figure(7); ax7 = axes(fig7); hold(ax7, 'on');grid on;title('Force to Platen'); legend();xlabel('Time (s)'); ylabel('Force (kN)');
fig8 = figure(8);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra');xlim([1 30]);subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra');xlim([0.1 5]);
set(fig8, 'WindowState', 'maximized');
color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; % Define colors for lines 1/3 and 2/4

%% Finding Response Spectre of Ground

f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);
[picos_ddx_ground , picos_x_ground] = ResponseSpectrum( t_vector , ddx_ref, f_vector , 1);

figure(fig8);
subplot(121)
grid on;
legend();
hold on
plot(f_vector, picos_ddx_ground(:, 1),'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Ground');% - Normal

subplot(122)
grid on;
legend();
hold on
plot(f_vector, picos_x_ground(:, 1),'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Ground ');%- Normal


%%

% Third plot
axes(ax3); % Activate the existing axes
x_T = lsim(G_xT_xref/s^2 ,  ddx_ref ,t_vector,'foh');
x_erro = x_T-x_ref;
mse = mean(x_erro.^2);
plot(t_vector,x_ref,"DisplayName","Reference")
plot(t_vector,x_T,"DisplayName","MSE="+string(mse))

% 5th plot
axes(ax5); % Activate the existing axes
plot(t_vector,x_erro,"DisplayName","Default")

% Fourth plot
axes(ax4); % Activate the existing axes
ddx_T = lsim(G_xT_xref, ddx_ref , t_vector ,'foh');
ddx_erro = ddx_T-ddx_ref;
mse = mean(ddx_erro.^2);
plot(t_vector,ddx_ref,"DisplayName","Reference")
plot(t_vector,ddx_T,"DisplayName","MSE="+string(mse))

% 6th plot
axes(ax6); % Activate the existing axes
plot(t_vector,ddx_erro,"DisplayName","Default")

% Second plot
axes(ax2); % Activate the existing axes
i_sv = lsim(G_c,   x_ref-x_T  , t_vector,'foh');
plot(t_vector,i_sv,"DisplayName","Default")

% 7th  plot 
axes(ax7); % Activate the existing axes
F_p_isv = lsim(G_Fp_isv,   i_sv  , t_vector,'foh');
plot(t_vector,F_p_isv/1e3,"DisplayName","Default") %/1e3 to display as kN



%% Finding Response Spectre for table

[picos_ddx_table , picos_x_table ] = ResponseSpectrum(  t_vector , ddx_T , f_vector, 1 );

figure(fig8);
subplot(121)
hold on
mse = mean((picos_ddx_table-picos_ddx_ground).^2);
plot(f_vector, picos_ddx_table(:, 1),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName',  sprintf('Platform - MSE= %.2e', mse(1)));

subplot(122)
hold on
mse = mean((picos_x_table-picos_x_ground).^2);
plot(f_vector, picos_x_table(:, 1),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName', sprintf('Platform - MSE= %.2e', mse(1)));



%% Driver update
driver = ddx_ref;
inv_G_xT_xref = inv(G_xT_xref );
% Get impulse response of G_inv (pad with more time for accuracy)
t = 0:0.001:10;
[imp_resp, t_imp] = impulse(inv_G_xT_xref, t);
% Convolve input with impulse response (this is time-domain filtering)
driver_correction =  conv(x_erro, imp_resp, 'same'); %lsim(inv_G_xT_xref,  ddx_erro , t_vector );
iteration_gain = 1;
driver = driver + driver_correction;

% % Define G_inv in frequency domain
% [mag, phase, w] = bode(G_inv);
% mag = squeeze(mag);
% phase = squeeze(phase);
% H = mag .* exp(1i*deg2rad(phase));  % Frequency response of G_inv
% 
% % FFT of input signal
% U = fft(u);
% 
% % Make sure frequency points match (may require interpolation!)
% % Multiply in frequency domain
% Y = H .* U;
% 
% % IFFT to get time-domain output
% y = real(ifft(Y));

%%

% Third plot
axes(ax3); % Activate the existing axes
x_T_driver = lsim(G_xT_xref/s^2 ,  driver ,t_vector,'foh');
erro = x_T_driver-x_ref;
mse = mean(erro.^2);
plot(t_vector,x_T_driver,"DisplayName","MSE Driver="+string(mse))

% 5th plot
axes(ax5); % Activate the existing axes
plot(t_vector,erro,"DisplayName","Driver")

% Fourth plot
axes(ax4); % Activate the existing axes
ddx_T_driver = lsim(G_xT_xref, driver , t_vector ,'foh');
ddx_erro = ddx_T_driver-ddx_ref;
mse = mean(ddx_erro.^2);
plot(t_vector,ddx_T_driver,"DisplayName","Driver")

% 6th plot
axes(ax6); % Activate the existing axes
plot(t_vector,ddx_erro,"DisplayName","Driver")

% Second plot
axes(ax2); % Activate the existing axes
i_sv = lsim(G_c,   driver-x_T_driver  , t_vector,'foh');
plot(t_vector,i_sv,"DisplayName","Driver")

% 7th  plot 
axes(ax7); % Activate the existing axes
F_p_isv = lsim(G_Fp_isv,   i_sv  , t_vector,'foh');
plot(t_vector,F_p_isv/1e3,"DisplayName","Default") %/1e3 to display as kN

%% Save all figures after plotting
% Folder path where you want to save the images
folderName = sprintf('Sim_Res/m_i=%.1f,f_1=%.1f, f_2=%.1f,zeta1=%.2f,zeta2=%.2f',mass*1e-3,f1,f2,zeta1,zeta2);

% Check if the folder already exists
if ~exist(folderName, 'dir')
    % Create the folder if it doesn't exist
    mkdir(folderName);
end

saveas(fig1,fullfile(folderName,'Bode_of_G_xT_xref.png'));
saveas(fig2,fullfile(folderName,'Input_to_Servo.png'));
saveas(fig3,fullfile(folderName,'Platen_Displacement.png'));
saveas(fig4,fullfile(folderName,'Platen_Acceleration.png'));
saveas(fig5,fullfile(folderName,'Platen_Displacement_Tracking_Error.png'));
saveas(fig6,fullfile(folderName,'Platen_Acceleration_Tracking_Error.png'));
saveas(fig7,fullfile(folderName,'Force_to_Platen.png'));
saveas(fig8,fullfile(folderName,'Response_Spectra.png'));
 
saveas(fig1,fullfile(folderName,'Bode_of_G_xT_xref.fig'));
saveas(fig2,fullfile(folderName,'Input_to_Servo.fig'));
saveas(fig3,fullfile(folderName,'Platen_Displacement.fig'));
saveas(fig4,fullfile(folderName,'Platen_Acceleration.fig'));
saveas(fig5,fullfile(folderName,'Platen_Displacement_Tracking_Error.fig'));
saveas(fig6,fullfile(folderName,'Platen_Acceleration_Tracking_Error.fig'));
saveas(fig7,fullfile(folderName,'Force_to_Platen.fig'));
saveas(fig8,fullfile(folderName,'Response_Spectra.fig'));
