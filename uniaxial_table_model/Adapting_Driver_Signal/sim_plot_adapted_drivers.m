clear;clc;close all;
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'

launch_Adapt =0; % Set to 1 to lauch Adapt.exe
return_on = 1; % Set to 1 for execution to stop before adapting drivers, or set to 0 if the adapted drivers have already been generated

% Load target
folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_Elcentro\';
target = 'elcentro.tgt'; 
addpath(folder);
LTF_to_TXT_then_load(target,'InputFolder', folder)
t_step = time_vector(2);

scale = 0.4;
if scale ~= 1
    x_tgt_T   = scale*x_tgt_T;
    x_tgt_L   = scale*x_tgt_L;
    ddx_tgt_T = scale*ddx_tgt_T;
    ddx_tgt_L = scale*ddx_tgt_L;
end

max_abs_x_tgt_T = max( abs(x_tgt_T ))
max_abs_x_tgt_L = max( abs(x_tgt_L ))

% figure;hold on; grid; legend;
% plot(time_vector , x_tgt_T)
% %plot(time_vector , x_tgt_L)
% s = tf('s');
% int_int__scaled_ddx_tgt_T = lsim(1/s^2 , ddx_tgt_T , time_vector , 'zoh');
% plot(time_vector , int_int__scaled_ddx_tgt_T)
% %
% % figure;hold on; grid; legend;
% % plot(time_vector , ddx_tgt_T)
% % s = tf('s');
% % int_int_ddx_tgt_T = lsim(1/s^2 , ddx_tgt_T , time_vector , 'zoh');
% % plot(time_vector , int_int_ddx_tgt_T)
% % 
% figure;hold on; grid; legend;
% plot(time_vector , ddx_tgt_T , '*-')
% dd_x_tgt_T= secondDerivativeTime(x_tgt_T,t_step);
% plot(time_vector , dd_x_tgt_T , '*-'  )
% % % dd_x_tgt_T_3= secondDerivativeTime3(x_tgt_T,t_step);
% % % plot(time_vector , dd_x_tgt_T_3 , '*-'  )
% % dd_x_tgt_T_7= secondDerivativeTime7(x_tgt_T,t_step);
% % plot(time_vector , dd_x_tgt_T_7 , '*-'  )

% figure;hold on; grid; legend;
% plot(time_vector , ddx_tgt_T , '*-')
% dd_x_tgt_T= secondDerivativeTime(x_tgt_T,t_step);
% plot(time_vector , dd_x_tgt_T , '*-'  )
% % plot(time_vector , x_tgt_T,'.-')
% xlim([0 t_step*11])

% sum = x_tgt_T(7) - 2* x_tgt_T(8)+ x_tgt_T(9)  
% ddx_Tstep7 = sum/t_step^2
% 
% ddx_diff=diff(x_tgt_T,2);
% plot(ddx_diff,'*-')

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
x_L_tuned = lsim(G_xT_xref_tuned ,  x_tgt_L , time_vector,'zoh');
ddx_L_tuned = secondDerivativeTime(x_L_tuned , t_step);

%% Optimal control
clsys = compute_Optimal_Controller( AA , BB , CC , DD);
[x_T_LQI,~, ~] = lsim(clsys, x_tgt_T, time_vector,'zoh'); % Simulate the closed-loop response using lsim:
ddx_T_LQI = secondDerivativeTime(x_T_LQI,t_step);
[x_L_LQI, ~, ~] = lsim(clsys, x_tgt_L, time_vector,'zoh'); % Simulate the closed-loop response using lsim:
ddx_L_LQI = secondDerivativeTime(x_L_LQI,t_step);

%% Lauch Adapt.exe % note the empty quotes "" are the window title placeholder
if launch_Adapt
    cmd = sprintf('start "" "%s"', fullfile('C:','Users','afons','OneDrive - Universidade de Lisboa','Controlo de Plataforma Sismica','uniaxial_table_model','Adapting_Driver_Signal','Adapt.exe.lnk'));
    system(cmd);
    fprintf("Launched Adapt.exe, continuing script...\n \n ");
end
fprintf("\n \n Go to Adapt.exe and generate driver 0 (Click 'Adapt Init' button) \n \n ")
if return_on
    return;
end   % execution stops here; lines below wonnt run

%% Simulation using updated driver 0
name = target(1 : end-4); %#ok<UNRCH>
LTF_to_TXT_then_load( [ name, '_0.DRV' ] ,'InputFolder',folder)
x_T_acq_0 = lsim(G_xT_xref ,  x_drv_T_0 , time_vector,'zoh');
ddx_T_acq_0 = secondDerivativeTime(x_T_acq_0 , t_step);
% writeTXT_then_LTF(time_vector,x_T_acq_0,ddx_T_acq_0,folder,[ name, '_0.ACQ.txt' ]);
x_L_acq_0 = lsim(G_xT_xref ,  x_drv_L_0 , time_vector,'zoh');
ddx_L_acq_0 = secondDerivativeTime(x_L_acq_0 , t_step);
writeTXT_then_LTF(time_vector,[x_T_acq_0,x_L_acq_0],[ddx_T_acq_0,ddx_L_acq_0],folder,[ name, '_0.ACQ.txt' ]);
fprintf("\n \n Go to Adapt.exe and generate driver 1 (Click 'Process' button)\n \n")
if return_on
    return;
end   % execution stops here; lines below wonnt run

%% Simulation using updated driver 1
LTF_to_TXT_then_load( [ name, '_1.DRV' ] ,'InputFolder',folder)
x_T_acq_1 = lsim(G_xT_xref ,  x_drv_T_1 , time_vector,'zoh');
ddx_T_acq_1 = secondDerivativeTime(x_T_acq_1 , t_step);
%writeTXT_then_LTF(time_vector,x_T_acq_1,ddx_T_acq_1,folder,[ name, '_1.ACQ.txt' ]);
x_L_acq_1 = lsim(G_xT_xref ,  x_drv_L_1 , time_vector,'zoh');
ddx_L_acq_1 = secondDerivativeTime(x_L_acq_1 , t_step);
writeTXT_then_LTF(time_vector,[x_T_acq_1,x_L_acq_1],[ddx_T_acq_1,ddx_L_acq_1],folder, [ name, '_1.ACQ.txt' ]); 
fprintf("\n \n Go to Adapt.exe and generate driver 2 (Click 'Next Iteration' and then 'Process' button) \n \n")
if return_on
    return;
end   % execution stops here; lines below wonnt run

%% Simulation using updated driver 2
LTF_to_TXT_then_load( [ name, '_2.DRV' ] ,'InputFolder',folder)
x_T_acq_2 = lsim(G_xT_xref ,  x_drv_T_2 , time_vector,'zoh');
ddx_T_acq_2 = secondDerivativeTime(x_T_acq_2 , t_step);
% writeTXT_then_LTF(time_vector,x_T_acq_2,ddx_T_acq_2,folder,[ name, '_2.ACQ.txt' ]);
x_L_acq_2 = lsim(G_xT_xref ,  x_drv_L_2 , time_vector,'zoh');
ddx_L_acq_2 = secondDerivativeTime(x_L_acq_2 , t_step);
writeTXT_then_LTF(time_vector,[x_T_acq_2,x_L_acq_2],[ddx_T_acq_2,ddx_L_acq_2],folder, [ name, '_2.ACQ.txt' ]); 

%% Create Figures - Transversal
close all;

% Response Spectra settings
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);

% Finding Response Spectre  of Target
[picos_ddx_tgt_T , picos_x_tgt_T] = ResponseSpectrum( time_vector , x_tgt_T , ddx_tgt_T, f_vector , 1);
[picos_ddx_tgt_L , picos_x_tgt_L] = ResponseSpectrum( time_vector , x_tgt_L , ddx_tgt_L, f_vector , 1);

% Response Spectre  of Optimal
[picos_ddx_T_tuned , picos_x_T_tuned] = ResponseSpectrum( time_vector , x_T_tuned , ddx_T_tuned, f_vector , 1);
[picos_ddx_L_tuned , picos_x_L_tuned] = ResponseSpectrum( time_vector , x_L_tuned , ddx_L_tuned, f_vector , 1);

% Response Spectre  of Optimal
[picos_ddx_T_LQI , picos_x_T_LQI] = ResponseSpectrum( time_vector , x_T_LQI , ddx_T_LQI, f_vector , 1);
[picos_ddx_L_LQI , picos_x_L_LQI] = ResponseSpectrum( time_vector , x_L_LQI , ddx_L_LQI, f_vector , 1);

% Computing Response spectra of Adapted
[picos_ddx_T_acq_0  , picos_x_T_acq_0 ] = ResponseSpectrum( time_vector , x_T_acq_0 , ddx_T_acq_0, f_vector , 1);
[picos_ddx_T_acq_1  , picos_x_T_acq_1 ] = ResponseSpectrum( time_vector , x_T_acq_1 , ddx_T_acq_1, f_vector , 1);
[picos_ddx_T_acq_2  , picos_x_T_acq_2 ] = ResponseSpectrum( time_vector , x_T_acq_2 , ddx_T_acq_2, f_vector , 1);
[picos_ddx_L_acq_0  , picos_x_L_acq_0 ] = ResponseSpectrum( time_vector , x_L_acq_0 , ddx_L_acq_0, f_vector , 1);
[picos_ddx_L_acq_1  , picos_x_L_acq_1 ] = ResponseSpectrum( time_vector , x_L_acq_1 , ddx_L_acq_1, f_vector , 1);
[picos_ddx_L_acq_2  , picos_x_L_acq_2 ] = ResponseSpectrum( time_vector , x_L_acq_2 , ddx_L_acq_2, f_vector , 1);

baseFolder = folder;   % Base folder where you want to create the timestamped subfolder
ts = datestr(now, 'yyyymmdd_HHMM');  % Create a timestamp string, e.g. '20250709_1530'
timeDir = fullfile(baseFolder, ts);  % Build the full path to the new folder
if ~exist(timeDir, 'dir')% Create it if it doesn't already exist
    mkdir(timeDir)
end

fig8 = figure(8);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Normal');xlim([1 20]);ylim([0 ceil(max(picos_ddx_T_tuned(1:385,1))) ])
subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Normal');xlim([0.1 5]);
color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; color4 = 'black';% Define colors for lines 1/3 and 2/4

figure(fig8); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_tgt_T,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
plot(f_vector, picos_ddx_T_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_T_tuned).^2 )));% - Normal
plot(f_vector, picos_ddx_T_LQI,'--', 'LineWidth' , 2 , 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e',   mean((picos_ddx_tgt_T-picos_ddx_T_LQI).^2 )));
plot(f_vector, picos_ddx_T_acq_0 ,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_0).^2 )));
plot(f_vector, picos_ddx_T_acq_1 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 1 -  MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_1).^2 )));
plot(f_vector, picos_ddx_T_acq_2 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 2 -  MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_2).^2 )));

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt_T,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
plot(f_vector, picos_x_T_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF - MSE= %.2e',     mean((picos_x_tgt_T-picos_x_T_tuned).^2 )));%- Normal
plot(f_vector, picos_x_T_LQI,'--', 'LineWidth' , 2, 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e',  mean((picos_x_tgt_T-picos_x_T_LQI).^2 )));
plot(f_vector, picos_x_T_acq_0, '-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_0).^2 )));
plot(f_vector, picos_x_T_acq_1, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 1 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_1).^2 )));
plot(f_vector, picos_x_T_acq_2, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 2 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_2).^2 )));

% Create Figures - Longitudinal
fig9 = figure(9);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Parallel');xlim([1 20]);ylim([0 ceil(max(picos_ddx_L_tuned(1:385,1))) ])
subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Parallel');xlim([0.1 5]);

figure(fig9); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_tgt_L,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
plot(f_vector, picos_ddx_L_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF -  MSE= %.2e',      mean((picos_ddx_tgt_L-picos_ddx_L_tuned).^2 )));% - Normal
plot(f_vector, picos_ddx_L_LQI,'--', 'LineWidth' , 2 , 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e',   mean((picos_ddx_tgt_L-picos_ddx_L_LQI).^2 )));
plot(f_vector, picos_ddx_L_acq_0 ,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_ddx_tgt_L-picos_ddx_L_acq_0).^2 )));
plot(f_vector, picos_ddx_L_acq_1 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 1 -  MSE= %.2e', mean((picos_ddx_tgt_L-picos_ddx_L_acq_1).^2 )));
plot(f_vector, picos_ddx_L_acq_2 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 2 -  MSE= %.2e', mean((picos_ddx_tgt_L-picos_ddx_L_acq_2).^2 )));

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt_L,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
plot(f_vector, picos_x_L_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF - MSE= %.2e',     mean((picos_x_tgt_L-picos_x_L_tuned).^2 )));%- Normal
plot(f_vector, picos_x_L_LQI,'--', 'LineWidth' , 2, 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e',  mean((picos_x_tgt_L-picos_x_L_LQI).^2 )));
plot(f_vector, picos_x_L_acq_0, '-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_x_tgt_L-picos_x_L_acq_0).^2 )));
plot(f_vector, picos_x_L_acq_1, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 1 - MSE= %.2e', mean((picos_x_tgt_L-picos_x_L_acq_1).^2 )));
plot(f_vector, picos_x_L_acq_2, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 2 - MSE= %.2e', mean((picos_x_tgt_L-picos_x_L_acq_2).^2 )));

%% Simulation using updated driver 3
% LTF_to_TXT_then_load( [ name, '_3.DRV' ] ,'InputFolder',folder)
% x_T_acq_3 = lsim(G_xT_xref ,  x_drv_T_3 , time_vector,'zoh');
% ddx_T_acq_3 = secondDerivativeTime(x_T_acq_3 , t_step);
% % writeTXT_then_LTF(time_vector,x_T_acq_3,ddx_T_acq_3,folder,[ name, '_3.ACQ.txt' ]);
% x_L_acq_3 = lsim(G_xT_xref ,  x_drv_L_3 , time_vector,'zoh');
% ddx_L_acq_3 = secondDerivativeTime(x_L_acq_3 , t_step);
% writeTXT_then_LTF(time_vector,[x_T_acq_3,x_L_acq_3],[ddx_T_acq_3,ddx_L_acq_3],folder, [ name, '_3.ACQ.txt' ]); 
% 
% [picos_ddx_T_acq_3  , picos_x_T_acq_3 ] = ResponseSpectrum( time_vector , x_T_acq_3 , ddx_T_acq_3, f_vector , 1);
% % Create Figures - Trnaversal
% figure(fig8); subplot(121);
% plot(f_vector, picos_ddx_T_acq_3 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 3 -  MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_3).^2 )));
% subplot(122);
% plot(f_vector, picos_x_T_acq_3, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 3 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_3).^2 )));
% 
% [picos_ddx_L_acq_3  , picos_x_L_acq_3 ] = ResponseSpectrum( time_vector , x_L_acq_3 , ddx_L_acq_3, f_vector , 1);
% % Create Figures - Longitudinal
% figure(fig9); subplot(121);
% plot(f_vector, picos_ddx_L_acq_3 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 3 -  MSE= %.2e', mean((picos_ddx_tgt_L-picos_ddx_L_acq_3).^2 )));
% subplot(122);
% plot(f_vector, picos_x_L_acq_3, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 3 - MSE= %.2e', mean((picos_x_tgt_L-picos_x_L_acq_3).^2 )));

%% Simulation using .acq from simulating directly in Adapt.exe
LTF_to_TXT_then_load( [ name, '_2.ACQ' ] , 'InputFolder', folder , 'OutputFolder', folder);

%%
set(fig8, 'WindowState', 'maximized');
exportgraphics(fig8,fullfile(timeDir,'Response_Spectra_N.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');

set(fig9, 'WindowState', 'maximized');
exportgraphics(fig9,fullfile(timeDir,'Response_Spectra_P.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');