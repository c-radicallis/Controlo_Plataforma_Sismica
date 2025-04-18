close all;clear;clc

% Loads variables and computes an optimal state-feedback control law for
% the tracking loop before running simulink
A = [0 1 0;
     0 0 1;
     1 0 0];    
B = [0.3; 1; -0.3];
C = [0 1 0];  
D = 0;
sys = ss(A, B, C, D);

nx = size(A,1);  
ny = size(C,1);

% --- Design the LQI Controller ---
% We design the LQI controller for the original plant.
Q = blkdiag(eye(nx), eye(ny))*1e3;
R = eye(size(B,2))*1e-3;  % Note: size(B,2) is the number of control inputs (should be 1)
K = lqi(ss(A, B, C, D), Q, R);
% K is the state-feedback gain that computes the control action.

% Define process noise matrices
G = eye(nx);         % Process noise matrix (assumed affecting all states)
H = zeros(ny, nx);   % No direct noise feedthrough
% Define noise covariance data.
% Qn: process noise covariance (nx-by-nx)
% Rn: measurement noise covariance (ny-by-ny)
Qn = eye(nx);
Rn = eye(ny);


%% --- Simulation ---
% Label plant inputs/outputs for interconnection:
sys.InputName = {'i_sv'};   % plant input: control signal
sys.OutputName = {'xT'};  % plant output
% --- Construct the Kalman Estimator ---

% Create the augmented system
sys_aug = ss(A, [B G], C, [D H]);
kest = kalman(sys_aug, Qn, Rn);

% --- Build the LQG Tracking Controller ---
% Combine the estimator and state-feedback gain into the tracking controller.
trksys = lqgtrack(kest, K);
% Label the controller’s I/O:
% The controller expects:
%   1st input: reference (r)
%   2nd input: measured output from the plant (y)
% And it produces:
%   output: control action (u)
trksys.InputName = {'r', 'xT'};
trksys.OutputName = {'i_sv'};
% --- Close the Loop ---
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
clsys = connect(sys, trksys, {'r'}, {'xT'})

% Define simulation time and reference signal (a step input of 1)
t = 0:0.01:15;          % time vector from 0 to 10 seconds
r = ones(length(t), 1); % step reference

% Simulate the closed-loop response using lsim:
[y_out, t_out, x] = lsim(clsys, r, t);

% Plot the response:
figure(1)
hold on
plot(t,r,'-.')
plot(t_out, y_out, 'LineWidth', 2)
plot(t_out,x(:,1:3))
plot(t_out,x(:,4:6),'--')
plot(t_out,x(:,7),':')
xlabel('Time (s)')
ylabel('Output y')
title('Closed-Loop Response with LQG Tracking Controller')
grid on

%  Load seismic signal and scale down if necessary
dados = load('elcentro.txt');
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
[y_out, t_out, x] = lsim(clsys, x_ref, t_vector);

% Plot the response:
figure(2)
hold on
plot(t_vector,x_ref,'-.')
plot(t_out, y_out, 'LineWidth', 2)
plot(t_out,x(:,1:3))
plot(t_out,x(:,4:6),'--')
plot(t_out,x(:,7),':')
xlabel('Time (s)')
ylabel('Output y')
title('Closed-Loop Response with LQG Tracking Controller')
grid on
