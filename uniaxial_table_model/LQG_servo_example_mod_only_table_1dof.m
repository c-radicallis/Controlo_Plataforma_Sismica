%clear;
close all;
clc;

% Structure parameters
mT=1.9751*1e3;       %Platen mass (mp=1.9751 t)
cT=5.78*1e3;        %Total damping, actuator + platen (ct=5.78 kN s/m1)

mass=2e3;

% 1st mode
m1 = mass; % kg
f1 = 0.2; % Hz   % 1.5 < f1 < 4
zeta1 = 0.1 ; % 2 < zeta1 < 10
%2nd mode
m2 = mass; % kg
f2 =10; % Hz % 6 < f2 < 10
zeta2 = 0.25; % 5 < zeta2 < 25

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
    (m1*k2 + m2*(k1+k2))/(2*m1*m2) - 0.5*sqrt( ((m1*k2+m2*(k1+k2))/(m1*m2))^2 - 4*( k1*k2 )/(m1*m2) ) == (2*pi*f1)^2,...
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
AA = [-1/tau_sv, 0              , 0 , 0      ;% Q_sv
             k_h/A , -k_h*k_pl/(A^2), 0 , -k_h   ;% Fp
                 0 ,              0 , 0 , 1      ;% xT
                 0 ,           1/mT , 0 , cT/mT];% dxT
          
BB = [k_svk_q/tau_sv ; zeros(3,1)];
CC = [0,0, 1 , 0];  % measuring xT
DD = 0;
sys = ss(AA,BB,CC,DD);
% Label plant inputs/outputs for interconnection:
sys.InputName = {'i_sv'};   % plant input: control signal
sys.OutputName = {'xT'};  % plant output

%
obs = obsv(AA, CC);
r_obsv = rank(obs)
ctrlb =ctrb(AA,BB);
r_ctrlb = rank(ctrlb)

format short g
obs=double(obs)
ctrlb = double(ctrlb)
format("default")

%% --- Build the Augmented System for Estimator Design ---
nx = size(AA,1);  
ny = size(CC,1);
% Define process noise matrices
G = eye(nx);         % Process noise matrix (assumed affecting all states)
H = zeros(ny, nx);   % No direct noise feedthrough

% Create the augmented system
sys_aug = ss(AA,[BB G], CC,[DD H]);
% Assign input groups:
%   Channel 1: control input (B)
%   Channels 2 to (nx+1): noise (G)
sys_aug.InputGroup.control = 1;
sys_aug.InputGroup.noise   = 2:(nx+1);
% Explicitly mark the control channel as KnownInput
sys_aug.InputGroup.KnownInput = 1;

%% --- Design the LQI Controller ---
% We design the LQI controller for the original plant.
Q = blkdiag(eye(nx), eye(ny));
R = 1e3*eye(size(BB,2));  % Note: size(B,2) is the number of control inputs (should be 1)
K = lqi(sys, Q, R)
% K is the state-feedback gain that computes the control action.

%% --- Construct the Kalman Estimator ---
% Define noise covariance data.
% Qn: process noise covariance (nx-by-nx)
% Rn: measurement noise covariance (ny-by-ny)
Qn = eye(nx);
Rn = eye(ny);
kest = kalman(sys_aug, Qn, Rn);

%% --- Build the LQG Tracking Controller ---
% Combine the estimator and state-feedback gain into the tracking controller.
trksys = lqgtrack(kest, K,'1dof');
% Label the controller’s I/O:
% The controller expects:
%   1st input: reference (r)
%   2nd input: measured output from the plant (y)
% And it produces:
%   output: control action (u)
trksys.InputName = {'e'};
trksys.OutputName = {'i_sv'};

Erro_soma = sumblk("e = x_ref - xT");

%% --- Close the Loop ---
% Now, interconnect the plant (sys) and the controller (trksys)
% using the connect function.
%
% The interconnection works as follows:
%  - The controller takes as input the reference 'r' and the plant output 'y'.
%  - The controller outputs 'u', which drives the plant.
%
% Define the connection:
%   External input: 'r'
%   External output: 'y'
clsys = connect(sys, trksys, Erro_soma, {'x_ref'}, {'xT'});

%% --- Simulation ---
% Define simulation time and reference signal (a step input of 1)
% t = 0:0.001:60;          % time vector from 0 to 10 seconds
% r = 10e-3*ones(length(t), 1); % step reference

%%  Load seismic signal and scale down if necessary
dados = load('elcentro.txt');
t_vector = dados(:,1);
t_step = t_vector(2);
ddx_ref = dados(:,2);

% Limits
lim_displacement = 100; % mm
lim_velocity = 0.4; % m/s
lim_force = 200e3; % N

% Scaling down if necessary
% displacements in milimeters
x_ref = lsim(1/s^2,  ddx_ref , t_vector ,'foh');
max_xref = max(x_ref);
scale=1;
while max_xref > lim_displacement
    scale = 0.95*scale;
    ddx_ref = 0.95*ddx_ref;
    x_ref = lsim(1e3/s^2,  ddx_ref , t_vector ,'foh');
    max_xref = max(x_ref);
end
scale
max_xref
v_ref =  lsim(1/s,  ddx_ref , t_vector ,'foh');
max_vref = max(v_ref)


% Simulate the closed-loop response using lsim:
[y_out, t_out, x] = lsim(clsys, x_ref, t_vector);

% Plot the response:
%figure
hold on
plot(t_vector,x_ref)
plot(t_out, y_out, 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Output y')
title('Closed-Loop Response with LQG Tracking Controller')
grid on
