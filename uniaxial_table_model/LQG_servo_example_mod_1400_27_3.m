close all;clear;clc

A = [0 1 0;
     0 0 1;
     1 0 0];
B = [0.3; 1; -0.3];
C = [0 1 0];   % measured output is the second state only
D = 0;
sys = ss(A, B, C, D)
% Label the plant’s I/O:
sys.InputName = {'i_sv'};   % plant input (control action)
sys.OutputName = {'xT'};  % plant output (only x2)

%% --- Build the Augmented System for Estimator Design ---
% Here, we use the same A and B, but the output remains C = [0 1 0].
nx = size(A,1);
ny = size(C,1);  % now ny = 1
% Define process noise matrices:
G = eye(nx);         % process noise input (affects all states)
H = zeros(ny, nx);   % no direct noise feedthrough (1-by-3 zero matrix)
sys_aug = ss(A, [B G], C, [D H]);
% Specify input groups:
%   Channel 1 is the control input,
%   Channels 2 to (nx+1) are noise.
sys_aug.InputGroup.control = 1;
sys_aug.InputGroup.noise   = 2:(nx+1);
% Indicate that only the first input is the known control input:
sys_aug.InputGroup.KnownInput = 1;

%% --- Design the LQI Controller ---
% Here we design an LQI controller based on the original plant.
% (For simplicity we use the same Q and R as before.)
Q = blkdiag(eye(nx), eye(ny));  % Note: this augments the state with the output error integrator.
R = eye(size(B,2));             % should be 1-by-1
K = lqi(ss(A, B, C, D), Q, R)

%% --- Design the Kalman Estimator ---
% Define noise covariances (dimensions must match sys_aug):
Qn = eye(nx);   % process noise covariance (3-by-3)
Rn = eye(ny);   % measurement noise covariance (1-by-1)
kest = kalman(sys_aug, Qn, Rn)
% (kest accepts two inputs: the known control input and the measured output.)

%% --- Form the LQG Tracking Controller ---
% Combine the estimator and state-feedback gain.
trksys = lqgtrack(kest, K)
% Specify its inputs:
%   First input: reference (r)
%   Second input: measured output (y)
trksys.InputName = {'r','xT'};
trksys.OutputName = {'i_sv'};
% Do not reassign the OutputName; we use the default.
% (Typically, the controller outputs its computed control action.)

%% --- Close the Loop ---
% Interconnect the plant and controller:
%  - The controller’s (only) output drives the plant.
%  - The plant output y (which is x2) is fed back to the controller.
% The free external input is the reference signal r.
%
% Use the plant output name ("y") and the controller’s default output names.
%ctrlOutNames = trksys.OutputName;  % e.g., {"OutputEstimate"} if only one output exists.
%clsys = connect(sys, trksys, {'r'}, [{'y'}, ctrlOutNames]);
clsys = connect(sys, trksys, {'r'}, {'xT'})

%% --- Simulate the Closed-Loop System ---
t = 0:0.01:10;           % time vector (0 to 10 seconds)
r = ones(length(t), 1);   % step reference (magnitude 1)

% Simulate the closed-loop response.
[Y, t_out] = lsim(clsys, r, t);
% The external outputs of clsys are:
%   Column 1: Plant output ("y")  – which is now x2 only
%   Column 2: Controller output (control action)

%% --- Simulate the Observer (Kalman Estimator) Separately ---
% The estimator (kest) accepts two inputs:
%   1. Known control input (u) – from the controller.
%   2. Measured plant output (y) – which is x2.
%
% Create the input matrix for the estimator.
% We assume the controller’s output is in column 2 of Y.
U_est = [Y(:,2), Y(:,1)];  % Order: [u, y]

% Simulate the estimator to obtain its state estimates:
[Ye, t_est, Xe] = lsim(kest, U_est, t);
% Xe contains the estimated states (all 3 states).

%% --- Plot the Results ---
figure

subplot(3,1,1)
plot(t_out, Y(:,1), 'b','LineWidth',2)
xlabel('Time (s)')
ylabel('y (Measured Output, x2)')
title('Plant Measured Output vs. Time')
grid on

subplot(3,1,2)
plot(t_out, Y(:,2), 'r','LineWidth',2)
xlabel('Time (s)')
ylabel('Control Action (u)')
title('Control Action vs. Time')
grid on

subplot(3,1,3)
plot(t_est, Xe(:,1), 'k--','LineWidth',1.5), hold on
plot(t_est, Xe(:,2), 'm--','LineWidth',1.5)
plot(t_est, Xe(:,3), 'g--','LineWidth',1.5)
xlabel('Time (s)')
ylabel('Observer States')
title('Observer (Kalman Filter) State Estimates vs. Time')
legend('x_1','x_2','x_3')
grid on
