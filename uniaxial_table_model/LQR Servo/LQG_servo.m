%clear;
clc;

mT=1.9751*1e3; %Platen mass (mp=1.9751 t)
cT=5.78*1e3;   %Total damping, actuator + platen (ct=5.78 kN s/m1)
mass=2e3;
m1 = mass; % kg % 1st mode
f1 = 1.5; % Hz   % 1.5 < f1 < 4
zeta1 = 0.02 ; % 2 < zeta1 < 10
m2 = mass; % kg %2nd mode
f2 =6; % Hz % 6 < f2 < 10
zeta2 = 0.05; % 5 < zeta2 < 25r

c1 = zeta1*2*m1*2*pi*f1; %N/m/s %coupled 2DOF system
c2 = zeta2*2*m2*2*pi*f2; %N/m/s
syms k1 k2
assume(k1 ,"positive")
assume(k2 ,"positive")
Y = vpasolve([
    (m1*k2 + m2*(k1+k2))/(2*m1*m2) - 0.5*sqrt( ((m1*k2+m2*(k1+k2))/(m1*m2))^2 - 4*( k1*k2 )/(m1*m2) ) == (2*pi*f1)^2  ,
    (m1*k2 + m2*(k1+k2))/(2*m1*m2)+ 0.5*sqrt( ((m1*k2+m2*(k1+k2))/(m1*m2))^2 - 4*( k1*k2 )/(m1*m2) ) == (2*pi*f2)^2,
], [k1,k2]);
k1 = double(Y.k1);
k2 = double(Y.k2);

tau_sv=0.0246     ; %Valve time constant (tsv=0.0246 s) % All converted to SI units
k_svk_q=1934.5*(1e-2)^3   ; %Valve flow gain (ksv∙kq=1934.5 cm3/s/V)
k_pl=1.67401e-7/1e3 ; %Valve pressure gain & leakadge factor (kpl=1.67401e-7 m3/s/kPa) %from Gidewon k_pl =  k_c + C_l 
Be=193716.28*1e3 ;  %Oil Bulk modulus (Be=193716.28 kPa)
Vt=0.002659   ;  %Oil Volume on actuator chamber (Vt=0.002659 m3)
A=0.012456    ;  %Piston area (A=0.012456 m2)
k_h=4*Be*A^2/Vt*1e3; %(kPa m1)

k_p=1.2993/1e-2; %SI units %Pgain (kp=1.2993 V/cm) 
G_c = tf(k_p,1);% Controller

s=tf('s'); % TF's
G_T=mT*s^2+(cT+c1)*s+k1;
G_1=m1*s^2+(c1+c2)*s+k1+k2;
G_2=m2*s^2+c2*s+k2;
G_T1=c1*s+k1;
G_21=c2*s+k2;
G_svq = k_svk_q/(1+tau_sv*s);
G_csv=G_c*G_svq;
G_x2_x1= G_21/G_2;
G_x1_xT = G_T1*G_2/(G_1*G_2-G_21^2);
G_xT_Fp = (G_1*G_2-G_21^2)/(G_T*G_1*G_2-G_T*G_21^2-G_2*G_T1^2);
G_Fp_xref = G_csv/( k_pl/A + A*s/k_h +G_xT_Fp*(G_csv + A*s));
G_xT_xref = G_Fp_xref * G_xT_Fp;
G_x1_xref = G_x1_xT*G_xT_xref;
G_x2_xT = G_x2_x1 * G_x1_xT;
G_Fp_isv = A*G_svq/( k_pl+A^2*s/k_h+A^2*s*G_xT_Fp );

AA = [-1/tau_sv, 0              , 0      , 0          , 0     , 0              , 0         , 0     ;
          k_h/A , -k_h*k_pl/(A^2), 0      , 0          , 0     , -k_h           , 0         , 0     ;
              0 ,              0 , 0      , 0          , 0     , 1              , 0         , 0     ;
              0 ,              0 , 0      , 0          , 0     , 0              , 1         , 0     ;
              0 ,              0 , 0      , 0          , 0     , 0              , 0         , 1     ;
              0 ,           1/mT , -k1/mT , k1/mT      ,  0    , (-cT - c1) /mT ,  c1 /mT   , 0     ;
              0 ,            0   , k1/m1  ,(-k1-k2)/m1 , k2/m1 , c1/m1          ,(-c1-c2)/m1, c2/m1 ;
              0 ,            0   , 0      , k2/m2      , -k2/m2, 0              , c2/m2     , -c2/m2];        
BB = [k_svk_q/tau_sv ; zeros(7,1)];
CC = [zeros(1,2), 1 , zeros(1,5)];  % measuring xT
DD = 0;
sys = ss(AA,BB,CC,DD);
sys.InputName = {'i_sv'};   % plant input: control signal
sys.OutputName = {'xT'};  % plant output

%% 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\
%load('mat and fig files\obsv_ctrb_vpa1e5.mat')

nx = size(AA,1);    % Number of states
nu = size(BB,2);    % Number of control inputs (should be 1)
ny = size(CC,1);    % Number of outputs

Q = 1e3*diag([zeros(1,nx),1]);%blkdiag(eye(nx), eye(ny));
R = 1e-9*eye(nu);
K = lqi(sys, Q, R)% Design the LQI controller for the original system

G = eye(nx);         % Process noise matrix (assuming full-state noise)
H = zeros(ny, nx);   % No direct noise feedthrough
sys_aug = ss(AA,[BB G], CC,[DD H]); % Create the augmented system
Qn = eye(nx);% Define noise covariance data
Rn = eye(ny);% Here Qn should be for process noise (nx-by-nx) and Rn for measurement noise (ny-by-ny)
kest = kalman(sys_aug, Qn, Rn)% Construct the Kalman estimator using the augmented system

trksys = lqgtrack(kest, K) % Build the LQG Tracking Controller- Combine the estimator and state-feedback gain into the tracking controller
trksys.InputName = {'x_ref', 'xT'};
trksys.OutputName = {'i_sv'};

clsys = connect(sys, trksys, {'x_ref'}, {'xT'}) %% --- Close the Loop ---

%% --- Simulation --
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
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

[x_T_LQG, t_out, x] = lsim(clsys, x_ref, t_vector,'foh'); % Simulate the closed-loop response using lsim:

figure(2); hold on; grid on; legend; % Plot the response:
plot(t_vector,x_ref,'-.')
plot(t_out, x_T_LQG, 'LineWidth', 2)
% plot(t_out,x(:,1:nx))
% plot(t_out,x(:,nx+1:end-1),'--')
% plot(t_out,x(:,end),':')
xlabel('Time (s)')
ylabel('Output y')
title('Closed-Loop Response with LQG Tracking Controller')

%% Plots
% clear; load('C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\LQR Servo\LQR_servo_1.mat'); x_T_LQG = y_out;

close all;
fig1 = figure(1); ax1 = axes(fig1); hold(ax1, 'on'); opts1=bodeoptions('cstprefs'); opts1.FreqUnits = 'Hz'; opts1.XLim={[1 40]}; % opts1.YLim={[-40 1]}; % opts1.MagVisible='off';
fig2 = figure(2); ax2 = axes(fig2); hold(ax2, 'on'); grid on; hold on; title('Input to Servo');  legend(); xlabel('Time (s)');  ylabel('Voltage (V)');
fig3 = figure(3); ax3 = axes(fig3); hold(ax3, 'on'); grid on; title(' Platen Displacement ');  hold on; legend(); xlabel('Time (s)');  ylabel('Displacement (mm)');
fig4 = figure(4); ax4 = axes(fig4); hold(ax4, 'on'); grid on; title('Platen Acceleration');  hold on; legend(); xlabel('Time (s)');  ylabel('Acceleration (m/s^2)'); 
fig5 = figure(5); ax5 = axes(fig5); hold(ax5, 'on');grid on;title(' Platen Displacement Tracking Error '); legend();xlabel('Time (s)'); ylabel('Error (mm)');
fig6 = figure(6); ax6 = axes(fig6); hold(ax6, 'on');grid on;title(' Platen Acceleration Tracking Error '); legend();xlabel('Time (s)'); ylabel('Error (m/s^2)');
fig7 = figure(7); ax7 = axes(fig7); hold(ax7, 'on');grid on;title('Force to Platen'); legend();xlabel('Time (s)'); ylabel('Force (kN)');
fig8 = figure(8);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra');xlim([1 30]);subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra');xlim([0.1 5]);
color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; % Define colors for lines 1/3 and 2/4

%% Finding Response Spectre of Ground
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);
[picos_ddx_ground , picos_x_ground] = ResponseSpectrum( t_vector , ddx_ref, f_vector , 1);

figure(fig8); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_ground(:, 1),'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Ground');% - Normal

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_ground(:, 1),'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Ground ');%- Normal

%%
axes(ax1); % First plot % Activate the existing axes
bodeplot(G_xT_xref,opts1);

axes(ax3);% Third plot % Activate the existing axes
x_T = lsim(G_xT_xref/s^2 ,  ddx_ref ,t_vector,'foh');
erro = x_T-x_ref;
mse = mean(erro.^2);
plot(t_vector,x_ref,"DisplayName","Reference")
plot(t_vector,x_T,"DisplayName","MSE="+string(mse))

axes(ax5);% 5th plot % Activate the existing axes
plot(t_vector,erro,"DisplayName","Default")

axes(ax4); % Activate the existing axes
ddx_T = lsim(G_xT_xref, ddx_ref , t_vector ,'foh');
erro = ddx_T-ddx_ref;
mse = mean(erro.^2);
plot(t_vector,ddx_ref,"DisplayName","Reference")
plot(t_vector,ddx_T,"DisplayName","MSE="+string(mse))

axes(ax6); % Activate the existing axes
plot(t_vector,erro,"DisplayName","Default")

axes(ax2); % Activate the existing axes
i_sv = lsim(G_c,  x_ref-x_T  , t_vector,'foh');
plot(t_vector,i_sv,"DisplayName","Default")

axes(ax7); % Activate the existing axes
F_p_isv = lsim(G_Fp_isv,   i_sv  , t_vector,'foh');
plot(t_vector,F_p_isv/1e3,"DisplayName","Default") %/1e3 to display as kN

%% Finding Response Spectre for table
[picos_ddx_table , picos_x_table ] = ResponseSpectrum(  t_vector , ddx_T , f_vector, 1 );

figure(fig8); subplot(121); hold on;
mse = mean((picos_ddx_table-picos_ddx_ground).^2);
plot(f_vector, picos_ddx_table(:, 1),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName',  sprintf('Platform - MSE= %.2e', mse(1)));

subplot(122);hold on;
mse = mean((picos_x_table-picos_x_ground).^2);
plot(f_vector, picos_x_table(:, 1),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName', sprintf('Platform - MSE= %.2e', mse(1)));

%%
axes(ax1); % Activate the existing axes
bodeplot(clsys,opts1);
legend( 'Default'   ,'LQG');
title('Bode of G\_xT\_xref'); 
grid on;

axes(ax3); % Activate the existing axes
erro = x_T_LQG-x_ref;
mse = mean(erro.^2);
plot(t_vector,x_T_LQG,"DisplayName","LQG MSE="+string(mse))

axes(ax5); % Activate the existing axes
plot(t_vector,erro,"DisplayName","LQG")

axes(ax4); % Activate the existing axes
ddx_T_LQG = lsim(clsys, ddx_ref ,t_vector,'foh');
erro = ddx_T_LQG-ddx_ref;
mse = mean(erro.^2);
plot(t_vector,ddx_T_LQG,"DisplayName","LQG MSE="+string(mse))

axes(ax6); hold on;% Activate the existing axes
plot(t_vector,erro,"DisplayName","LQG")
%%
axes(ax2); hold on;
i_sv = lsim(trksys ,   [x_ref , x_T_LQG]  ,t_vector,'foh');
plot(t_vector,i_sv,"DisplayName","LQG")

axes(ax7); hold on;
F_p_isv = x(:,2); % 2nd element of state vector
plot(t_vector,F_p_isv*1e-3,"DisplayName","LQG")

%% Finding Response Spectre for table LQG
[picos_ddx_table_LQG , picos_x_table_LQG ] = ResponseSpectrum( t_vector , ddx_T_LQG, f_vector , 1 );

figure(fig8); subplot(121); hold on;
mse = mean((picos_ddx_table_LQG-picos_ddx_ground).^2);
plot(f_vector, picos_ddx_table_LQG(:, 1),'-', 'LineWidth' , 2, 'Color', color3, 'DisplayName', sprintf('LQG Platform - MSE= %.2e', mse(1)));

subplot(122);hold on;
mse = mean((picos_x_table_LQG-picos_x_ground).^2);
plot(f_vector, picos_x_table_LQG(:, 1),'-', 'LineWidth' , 2, 'Color', color3, 'DisplayName',  sprintf('LQG Platform - MSE= %.2e', mse(1)));

%% Save all figures after plotting
folderName = sprintf('Sim_Res_LQG/m_i=%.1f,f_1=%.1f, f_2=%.1f,zeta1=%.2f,zeta2=%.2f',mass*1e-3,f1,f2,zeta1,zeta2); % Folder path where you want to save the images
if ~exist(folderName, 'dir')% Check if the folder already exists
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

