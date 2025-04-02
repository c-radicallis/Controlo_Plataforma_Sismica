clear; clc

% Structure parameters

mT=1.9751*1e3;       %Platen mass (mp=1.9751 t)
cT=5.78*1e3;        %Total damping, actuator + platen (ct=5.78 kN s/m1)

mass=2e3;

% 1st mode
m1 = mass; % kg
f1 = 2; % Hz   % 1.5 < f1 < 4
zeta1 = 0.1 ; % 2 < zeta1 < 10
%2nd mode
m2 = mass; % kg
f2 =10; % Hz % 6 < f2 < 10
zeta2 = 0.06; % 5 < zeta2 < 25r

% Controller
k_p=1.2993/1e-2; %SI units %Pgain (kp=1.2993 V/cm) 
G_c = tf(k_p,1);

%coupled 2DOF system
c1 = zeta1*2*m1*2*pi*f1; %N/m/s
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

%Servo-valve parameters
% All converted to SI units
tau_sv=0.0246     ; %Valve time constant (tsv=0.0246 s)
k_svk_q=1934.5*(1e-2)^3   ; %Valve flow gain (ksv∙kq=1934.5 cm3/s/V)

%from Gidewon k_pl =  k_c + C_l 
k_pl=1.67401e-7/1e3 ; %Valve pressure gain & leakadge factor (kpl=1.67401e-7 m3/s/kPa)

Be=193716.28*1e3 ;  %Oil Bulk modulus (Be=193716.28 kPa)
Vt=0.002659   ;  %Oil Volume on actuator chamber (Vt=0.002659 m3)
A=0.012456    ;  %Piston area (A=0.012456 m2)
k_h=4*Be*A^2/Vt*1e3; %(kPa m1)

% mT=1.9751*1e3;       %Platen mass (mp=1.9751 t)
% cT=5.78*1e3;        %Total damping, actuator + platen (ct=5.78 kN s/m1)


% TF's
s=tf('s')  ;

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

% state space model
% digits(1e2)
% AA = vpa([-1/tau_sv, 0              , 0      , 0          , 0     , 0              , 0         , 0     ;
%           k_h/A , -k_h*k_pl/(A^2), 0      , 0          , 0     , -k_h           , 0         , 0     ;
%               0 ,              0 , 0      , 0          , 0     , 1              , 0         , 0     ;
%               0 ,              0 , 0      , 0          , 0     , 0              , 1         , 0     ;
%               0 ,              0 , 0      , 0          , 0     , 0              , 0         , 1     ;
%               0 ,           1/mT , -k1/mT , k1/mT      ,  0    , (-cT - c1) /mT ,  c1 /mT   , 0     ;
%               0 ,            0   , k1/m1  ,(-k1-k2)/m1 , k2/m1 , c1/m1          ,(-c1-c2)/m1, c2/m1 ;
%               0 ,            0   , 0      , k2/m2      , -k2/m2, 0              , c2/m2     , -c2/m2]);
% 
% BB = vpa([k_svk_q/tau_sv ; zeros(7,1)]);
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
% %%
% obs = vpa(obsv(AA, CC));
% r_obsv = rank(obs)
% ctrlb = vpa(ctrb(AA,BB));
% r_ctrlb = rank(ctrlb)
% 
% %%
% format short g
% obs=double(obs)
% ctrlb = double(ctrlb)
% format

%% 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\
% clear
% load('mat and fig files\obsv_ctrb_vpa1e5.mat')

% ss_model =
% 
%   A = 
%                x1          x2          x3          x4          x5          x6          x7          x8
%    x1      -40.65           0           0           0           0           0           0           0
%    x2    3.63e+12  -4.878e+04           0           0           0  -4.521e+10           0           0
%    x3           0           0           0           0           0           1           0           0
%    x4           0           0           0           0           0           0           1           0
%    x5           0           0           0           0           0           0           0           1
%    x6           0   0.0005063      -13.03       13.03           0      -3.435       0.509           0
%    x7           0           0       12.87      -187.2       174.4      0.5027      -2.765       2.262
%    x8           0           0           0       174.4      -174.4           0       2.262      -2.262
% 
%   B = 
%             u1
%    x1  0.07864
%    x2        0
%    x3        0
%    x4        0
%    x5        0
%    x6        0
%    x7        0
%    x8        0
% 
%   C = 
%        x1  x2  x3  x4  x5  x6  x7  x8
%    y1   0   0   1   0   0   0   0   0
% 
%   D = 
%        u1
%    y1   0

%%
% obs = vpa(obsv(AA, CC));
% r_obsv = rank(obs)
% obs=double(obs)
% ctrlb = vpa(ctrb(AA,BB));
% r_ctrlb = rank(ctrlb)
% ctrlb = double(ctrlb)

% r_obsv =
% 
%      8
% 
% 
% obs =
% 
%    1.0e+28 *
% 
%          0         0    0.0000         0         0         0         0         0
%          0         0         0         0         0    0.0000         0         0
%          0    0.0000   -0.0000    0.0000         0   -0.0000    0.0000         0
%     0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0000    0.0000
%    -0.0000    0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000
%     0.0000   -0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.0000
%    -0.0000    0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000
%     1.0118   -0.0000   -0.0000    0.0000   -0.0000   -0.0126    0.0000   -0.0000
% 
% 
% r_ctrlb =
% 
%      8
% 
%  ctrlb =
% 
% 1.0e+39 *
% 
% 0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000
%      0    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0001    3.6671
%      0         0         0    0.0000   -0.0000    0.0000   -0.0000    0.0000
%      0         0         0         0    0.0000   -0.0000    0.0000   -0.0000
%      0         0         0         0         0    0.0000   -0.0000    0.0000
%      0         0    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000
%      0         0         0    0.0000   -0.0000    0.0000   -0.0000    0.0000
%      0         0         0         0    0.0000   -0.0000    0.0000   -0.0000


%% 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\
%load('mat and fig files\obsv_ctrb_vpa1e5.mat')

nx = size(AA,1);    % Number of states
nu = size(BB,2);    % Number of control inputs (should be 1)
ny = size(CC,1);    % Number of outputs

% Define process noise matrices
G = eye(nx);         % Process noise matrix (assuming full-state noise)
H = zeros(ny, nx);   % No direct noise feedthrough

% Create the augmented system
sys_aug = ss(AA,[BB G], CC,[DD H]); 

% Assign input groups:
% % - Channel 1 is control,
% % - Channels 2 to (nx+1) are noise.
% sys_aug.InputGroup.control = 1;
% sys_aug.InputGroup.noise   = 2:(nx+1);
% % Explicitly set the KnownInput group to only the control input
% sys_aug.InputGroup.KnownInput = 1;

% Design the LQI controller for the original system
Q = 1e3*diag([zeros(1,nx),1]);%blkdiag(eye(nx), eye(ny));
R = 1e-9*eye(nu);
K = lqi(sys, Q, R)

% Define noise covariance data
% Here Qn should be for process noise (nx-by-nx) and Rn for measurement noise (ny-by-ny)
Qn = eye(nx);
Rn = eye(ny);

% Construct the Kalman estimator using the augmented system
kest = kalman(sys_aug, Qn, Rn)

trksys = lqgtrack(kest, K) % Build the LQG Tracking Controller --- Combine the estimator and state-feedback gain into the tracking controller
trksys.InputName = {'x_ref', 'xT'};
trksys.OutputName = {'i_sv'};

clsys = connect(sys, trksys, {'x_ref'}, {'xT'}) %% --- Close the Loop ---

%% --- Simulation ---

% % Define simulation time and reference signal (a step input of 1)
% t = 0:0.005:3;          % time vector from 0 to 10 seconds
% r = 0.1*ones(length(t), 1); % step reference
% 
% Simulate the closed-loop response using lsim:
% [y_out, t_out, x] = lsim(clsys, r, t);
% 
% Plot the response:
% figure(1)
% hold on
% plot(t,r,'-.')
% plot(t_out, y_out, 'LineWidth', 2)
% plot(t_out,x(:,1:nx))
% plot(t_out,x(:,nx+1:end-1),'--')
% plot(t_out,x(:,end),':')
% xlabel('Time (s)')
% ylabel('Output y')
% title('Closed-Loop Response with LQG Tracking Controller')
% grid on
% legend

%  Load seismic signal and scale down if necessary
dados = load('uniaxial_table_model\elcentro.txt');
t_vector = dados(:,1);
t_step = t_vector(2);
ddx_ref = dados(:,2);

% Limits
lim_displacement = 0.1; % m
lim_velocity = 0.4; % m/s
lim_force = 200e3; % N

s=tf('s')  ;
x_ref = lsim(1/s^2,  ddx_ref , t_vector ,'foh');
max_xref = max(x_ref);
scale=1;
% Scaling down if necessary
while max_xref > lim_displacement
    scale = 0.95*scale;
    ddx_ref = 0.95*ddx_ref;
    x_ref = lsim(1/s^2,  ddx_ref , t_vector ,'foh');
    max_xref = max(x_ref);
end
scale;
max_xref;
v_ref =  lsim(1/s,  ddx_ref , t_vector ,'foh');
max_vref = max(v_ref);


% Simulate the closed-loop response using lsim:
[x_T_LQG, t_out, x] = lsim(clsys, x_ref, t_vector,'foh');

% Plot the response:
figure(2)
hold on
plot(t_vector,x_ref,'-.')
plot(t_out, x_T_LQG, 'LineWidth', 2)
% plot(t_out,x(:,1:nx))
% plot(t_out,x(:,nx+1:end-1),'--')
% plot(t_out,x(:,end),':')
xlabel('Time (s)')
ylabel('Output y')
title('Closed-Loop Response with LQG Tracking Controller')
grid on
legend

%% Plots
% clear;
% load('C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\LQR Servo\LQR_servo_1.mat')
% x_T_LQG = y_out;

close all;

fig1 = figure(1); %clf; % Clear figure if needed
ax1 = axes(fig1); hold(ax1, 'on');
opts1=bodeoptions('cstprefs');
opts1.FreqUnits = 'Hz';
opts1.XLim={[1 40]};
% opts1.YLim={[-40 1]};
% opts1.MagVisible='off';

fig2 = figure(2); %clf;
ax2 = axes(fig2); hold(ax2, 'on');
grid on
hold on
title('Input to Servo'); 
legend();
xlabel('Time (s)'); 
ylabel('Voltage (V)');

fig3 = figure(3); %clf;
ax3 = axes(fig3); hold(ax3, 'on');
grid on
title(' Platen Displacement '); 
hold on
legend()
xlabel('Time (s)'); 
ylabel('Displacement (mm)');

fig4 = figure(4); %clf;
ax4 = axes(fig4); hold(ax4, 'on');
grid on
title('Platen Acceleration'); 
hold on
legend()
xlabel('Time (s)'); 
ylabel('Acceleration (m/s^2)');

fig5= figure(5); %clf; 
ax5 = axes(fig5); hold(ax5, 'on');
grid on
title(' Platen Displacement Tracking Error '); 
legend()
xlabel('Time (s)'); 
ylabel('Error (mm)');

fig6= figure(6); %clf; 
ax6 = axes(fig6); hold(ax6, 'on');
grid on
title(' Platen Acceleration Tracking Error '); 
legend()
xlabel('Time (s)'); 
ylabel('Error (m/s^2)');

fig7= figure(7); %clf; 
ax7 = axes(fig7); hold(ax7, 'on');
grid on
title('Force to Platen'); 
legend();
xlabel('Time (s)'); 
ylabel('Force (kN)');

fig8 = figure(8);%
subplot(121);
grid on;
xlabel('Frequency (Hz)');
ylabel('Acceleration (m/s^2)');
title('Acceleration Response Spectra');
xlim([1 30]);
subplot(122);
grid on;
xlabel('Frequency (Hz)');
ylabel('Displacement (m)');
title('Displacement Response Spectra');
xlim([0.1 5]);
% Define colors for lines 1/3 and 2/4
color1 = 'blue';
color2 = 'red' ;
color3 = '#EDB120';

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
% First plot
axes(ax1); % Activate the existing axes
bodeplot(G_xT_xref,opts1);

% Third plot
axes(ax3); % Activate the existing axes

x_T = lsim(G_xT_xref/s^2 ,  ddx_ref ,t_vector,'foh');
erro = x_T-x_ref;
mse = mean(erro.^2);
plot(t_vector,x_ref,"DisplayName","Reference")
plot(t_vector,x_T,"DisplayName","MSE="+string(mse))

%%
% 5th plot
axes(ax5); % Activate the existing axes
plot(t_vector,erro,"DisplayName","Default")

% Fourth plot
axes(ax4); % Activate the existing axes
ddx_T = lsim(G_xT_xref, ddx_ref , t_vector ,'foh');
erro = ddx_T-ddx_ref;
mse = mean(erro.^2);
plot(t_vector,ddx_ref,"DisplayName","Reference")
plot(t_vector,ddx_T,"DisplayName","MSE="+string(mse))

% 6th plot
axes(ax6); % Activate the existing axes
plot(t_vector,erro,"DisplayName","Default")

% Second plot
axes(ax2); % Activate the existing axes
i_sv = lsim(G_c,  x_ref-x_T  , t_vector,'foh');
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


%%
% First plot
axes(ax1); % Activate the existing axes
bodeplot(clsys,opts1);
legend( 'Default'   ,'LQG');
title('Bode of G\_xT\_xref'); 
grid on;

% displacements in milimeters
axes(ax3); % Activate the existing axes
erro = x_T_LQG-x_ref;
mse = mean(erro.^2);
plot(t_vector,x_T_LQG,"DisplayName","LQG MSE="+string(mse))

% 5th plot
axes(ax5); % Activate the existing axes
plot(t_vector,erro,"DisplayName","Tuned")

% Fourth plot
axes(ax4); % Activate the existing axes
ddx_T_LQG = lsim(clsys, ddx_ref ,t_vector,'foh');
erro = ddx_T_LQG-ddx_ref;
mse = mean(erro.^2);
plot(t_vector,ddx_T_LQG,"DisplayName","Tuned MSE="+string(mse))

% 6th plot
axes(ax6); % Activate the existing axes
hold on
plot(t_vector,erro,"DisplayName","Tuned")


% % Second plot
% axes(ax2); % Activate the existing axes
% hold on
% i_sv = lsim(G_c ,   x_ref-x_T_LQG  ,t_vector,'foh');
% plot(t_vector,i_sv,"DisplayName","Tuned")
%  plot
axes(ax7); % Activate the existing axes
hold on
F_p_isv = x(:,2);
plot(t_vector,F_p_isv*1e-3,"DisplayName","Tuned")



%% Finding Response Spectre for table tuned

[picos_ddx_table_tuned , picos_x_table_tuned ] = ResponseSpectrum( t_vector , ddx_T_LQG, f_vector , 1 );


figure(fig8);
subplot(121)
hold on
mse = mean((picos_ddx_table_tuned-picos_ddx_ground).^2);
plot(f_vector, picos_ddx_table_tuned(:, 1),'-', 'LineWidth' , 2, 'Color', color3, 'DisplayName', sprintf('Tuned Platform - MSE= %.2e', mse(1)));

subplot(122)
hold on
mse = mean((picos_x_table_tuned-picos_x_ground).^2);
plot(f_vector, picos_x_table_tuned(:, 1),'-', 'LineWidth' , 2, 'Color', color3, 'DisplayName',  sprintf('Tuned Platform - MSE= %.2e', mse(1)));


%% Save all figures after plotting
% Folder path where you want to save the images
folderName = sprintf('Sim_Res_LQG/m_i=%.1f,f_1=%.1f, f_2=%.1f,zeta1=%.2f,zeta2=%.2f',mass*1e-3,f1,f2,zeta1,zeta2);

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

