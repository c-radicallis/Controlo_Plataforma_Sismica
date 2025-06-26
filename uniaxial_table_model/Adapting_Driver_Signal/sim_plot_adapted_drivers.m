clear;clc;close all;
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'\Adapting_Driver_Signal\PRJ_project_1\

return_on = 0; % Set to 1 for execution to stop before adapting drivers, or set to 0 if the adapted drivers have already been generated

%% Load target
folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_project_1';
target = 'LAquilaReducedScale.tgt';   % or get from user input % 2. Define only the name (no folder); you can prompt the u
LTF_to_TXT_then_load(target,'InputFolder', folder)
t_step = time_vector(2);

%% Response Spectra settings
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);

%% Finding Response Spectre  of Target
[picos_ddx_tgt , picos_x_tgt] = ResponseSpectrum( time_vector , x_tgt_T , ddx_tgt_T, f_vector , 1);

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
[~,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ , AA , BB , CC , DD  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);


%% Finding Response Spectre of Table with tuned PID

tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
cutoff_frequency = 20; % Hz
G_c   = pidtune(G_Fp_isv*G_xT_Fp,'PIDF',cutoff_frequency*2*pi,tuner_opts)
[s,~,~,~,~ ,~ ,~,~,~,~,~,~,G_xT_xref_tuned,~,~ , ~ ,~,~,~,~ ,~  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);

x_T_tuned = lsim(G_xT_xref_tuned ,  x_tgt_T , time_vector,'zoh');
ddx_T_tuned = secondDerivativeTime(x_T_tuned , t_step);
[picos_ddx_tuned , picos_x_tuned] = ResponseSpectrum( time_vector , x_T_tuned , ddx_T_tuned, f_vector , 1);

%% Optimal control
sys = ss(AA,BB,CC,DD);
nx = size(AA,1);    % Number of states
nu = size(BB,2);    % Number of control inputs (should be 1)
ny = size(CC,1);    % Number of outputs

plant_aug = ss(AA, BB,[eye(nx);CC],DD);
plant_aug.InputName = {'i_sv'};   % plant input: control signal
plant_aug.OutputName = {'Qsv' , 'Fp' , 'xT' , 'x1' , 'x2','dxT' , 'dx1' , 'dx2' , 'y_xT'};  % plant output
sumblk1 = sumblk('e = x_tgt - y_xT'); % Compute the error signal: e = r - y
integrator = tf(1,[1 0]); % The integrator integrates the tracking error.
integrator.InputName = {'e'};    % error: e = r - y
integrator.OutputName = {'xi'};  % integrated error

Q = 1e3*diag([zeros(1,nx),1]);%blkdiag(eye(nx), eye(ny));
R = 1e-9*eye(nu);
K_lqi = lqi(sys, Q, R)% Design the LQI controller for the original system
K  = K_lqi(1:nx);      % state feedback gains
Ki = K_lqi(end);        % integrator gain
controller = ss([], [], [], -[K, Ki]); %   u = -[K  Ki] * [x; xi]
controller.InputName = {'Qsv' , 'Fp' , 'xT' , 'x1' , 'x2','dxT' , 'dx1' , 'dx2' , 'xi'};
controller.OutputName = {'i_sv'};
clsys = connect(plant_aug,  controller , integrator, sumblk1, 'x_tgt', 'y_xT')

[x_T_LQI, t_out, x] = lsim(clsys, x_tgt_T, time_vector,'zoh'); % Simulate the closed-loop response using lsim:
ddx_T_LQI = secondDerivativeTime(x_T_LQI,t_step);
[picos_ddx_LQI , picos_x_LQI] = ResponseSpectrum( time_vector , x_T_LQI , ddx_T_LQI, f_vector , 1);

%% Simulation using updated driver 0

LTF_to_TXT_then_load('LAquilaReducedScale_0.DRV')
x_acq_0 = lsim(G_xT_xref ,  x_drv_T_0 , time_vector,'zoh');
ddx_acq_0 = secondDerivativeTime(x_acq_0 , t_step);
writeTXT_then_LTF(time_vector,x_acq_0,ddx_acq_0,folder, 'LAquilaReducedScale_0.ACQ');
[picos_ddx_acq_0  , picos_x_acq_0 ] = ResponseSpectrum( time_vector , x_acq_0 , ddx_acq_0, f_vector , 1);

if return_on
    return;
end   % execution stops here; lines below wonnt run
%% Simulation using updated driver 1

LTF_to_TXT_then_load('LAquilaReducedScale_1.DRV')
x_acq_1 = lsim(G_xT_xref ,  x_drv_T_1 , time_vector,'zoh');
ddx_acq_1 = secondDerivativeTime(x_acq_1 , t_step);
[picos_ddx_acq_1  , picos_x_acq_1 ] = ResponseSpectrum( time_vector , x_acq_1 , ddx_acq_1, f_vector , 1);
writeTXT_then_LTF(time_vector,x_acq_1,ddx_acq_1,folder, 'LAquilaReducedScale_1.ACQ');

if return_on
    return;
end   % execution stops here; lines below won’t run
%% Simulation using updated driver 2
    
LTF_to_TXT_then_load('LAquilaReducedScale_2.DRV')
x_acq_2 = lsim(G_xT_xref ,  x_drv_T_2 , time_vector,'zoh');
ddx_acq_2 = secondDerivativeTime(x_acq_2 , t_step);
[picos_ddx_acq_2  , picos_x_acq_2 ] = ResponseSpectrum( time_vector , x_acq_2 , ddx_acq_2, f_vector , 1);
writeTXT_then_LTF(time_vector,x_acq_2,ddx_acq_2,folder, 'LAquilaReducedScale_2.ACQ');

%% Create Figures
fig8 = figure(8);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra');xlim([1 20]);subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra');xlim([0.1 5]);
color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; color4 = 'black';% Define colors for lines 1/3 and 2/4

figure(fig8); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_tgt,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
plot(f_vector, picos_ddx_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF -  MSE= %.2e', mean((picos_ddx_tgt-picos_ddx_tuned).^2 )));% - Normal
plot(f_vector, picos_ddx_LQI,'--', 'LineWidth' , 2 , 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e', mean((picos_ddx_tgt-picos_ddx_LQI).^2 )));
plot(f_vector, picos_ddx_acq_0 ,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_ddx_tgt-picos_ddx_acq_0).^2 )));
plot(f_vector, picos_ddx_acq_1 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 1 -  MSE= %.2e', mean((picos_ddx_tgt-picos_ddx_acq_1).^2 )));
plot(f_vector, picos_ddx_acq_2 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 2 -  MSE= %.2e', mean((picos_ddx_tgt-picos_ddx_acq_2).^2 )));

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
plot(f_vector, picos_x_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF - MSE= %.2e', mean((picos_x_tgt-picos_x_tuned).^2 )));%- Normal
plot(f_vector, picos_x_LQI,'--', 'LineWidth' , 2, 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e', mean((picos_x_tgt-picos_x_LQI).^2 )));
plot(f_vector, picos_x_acq_0, '-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver - MSE= %.2e', mean((picos_x_tgt-picos_x_acq_0).^2 )));
plot(f_vector, picos_x_acq_1, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 1 - MSE= %.2e', mean((picos_x_tgt-picos_x_acq_1).^2 )));
plot(f_vector, picos_x_acq_2, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 2 - MSE= %.2e', mean((picos_x_tgt-picos_x_acq_2).^2 )));

set(fig8, 'WindowState', 'maximized');
exportgraphics(fig8,fullfile('C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal','Response_Spectra.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');