clear;clc

%% --- Define the Original System (Plant) ---
A = [0 1 0;
     0 0 1;
     1 0 0];    
B = [0.3; 1; -0.3];
C = [0 1 0];  
D = 0;
sys = ss(A, B, C, D);
% Label plant inputs/outputs for interconnection:
sys.InputName = {'u'};   % plant input: control signal
sys.OutputName = {'y'};  % plant output

observability = rank(obsv(A,C))
controlability = rank(ctrb(A,B))

% digits(1e5)
% A = vpa([0 1 0;
%      0 0 1;
%      1 0 0]);    
% B = vpa([0.3; 0; -0.3]);
% C = vpa([0 1 0]);  
% D = 0;
% observability_vpa= rank(vpa(obsv(A,C)))
% controlability_vpa = rank(vpa(ctrb(A,B)))

%% --- Build the Augmented System for Estimator Design ---
nx = size(A,1);  
ny = size(C,1);
% Define process noise matrices
G = eye(nx);         % Process noise matrix (assumed affecting all states)
H = zeros(ny, nx);   % No direct noise feedthrough

% Create the augmented system
sys_aug = ss(A, [B G], C, [D H]);
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
R = eye(size(B,2));  % Note: size(B,2) is the number of control inputs (should be 1)
K = lqi(ss(A, B, C, D), Q, R);
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
trksys = lqgtrack(kest, K);
% Label the controller’s I/O:
% The controller expects:
%   1st input: reference (r)
%   2nd input: measured output from the plant (y)
% And it produces:
%   output: control action (u)
trksys.InputName = {'r', 'y'};
trksys.OutputName = {'u'};

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
%
clsys = connect(sys, trksys, {'r'}, {'y'});

%% --- Simulation ---
% Define simulation time and reference signal (a step input of 1)
t = 0:0.01:10;          % time vector from 0 to 10 seconds
r = ones(length(t), 1); % step reference

% Simulate the closed-loop response using lsim:
[y_out, t_out, x] = lsim(clsys, r, t);

% Plot the response:
figure
plot(t_out, y_out, 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Output y')
title('Closed-Loop Response with LQG Tracking Controller')
grid on
