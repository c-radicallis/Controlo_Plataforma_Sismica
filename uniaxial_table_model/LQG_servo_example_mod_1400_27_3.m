clear; clc

%% --- Define the Original System (Plant) ---
A = [0 1 0;
     0 0 1;
     1 0 0];    
B = [0.3; 1; -0.3];
C = [1.9 1.3 1];  
D = 0;
sys = ss(A, B, C, D);
% Label plant I/O for interconnection:
sys.InputName = {'u'};   % plant input (control action)
sys.OutputName = {'y'};  % plant output

%% --- Build the Augmented System for Estimator Design ---
nx = size(A,1);
ny = size(C,1);
% Define process noise matrices:
G = eye(nx);         % process noise matrix (affects all states)
H = zeros(ny, nx);   % no direct noise feedthrough
sys_aug = ss(A, [B G], C, [D H]);
% Specify input groups:
%   Channel 1 is the control input, channels 2 to (nx+1) are noise.
sys_aug.InputGroup.control = 1;
sys_aug.InputGroup.noise   = 2:(nx+1);
% Indicate that only the first input is a known control input:
sys_aug.InputGroup.KnownInput = 1;

%% --- Design the LQI Controller ---
Q = blkdiag(eye(nx), eye(ny));
R = eye(size(B,2));  % should be 1-by-1
K = lqi(ss(A, B, C, D), Q, R);

%% --- Design the Kalman Estimator ---
Qn = eye(nx);   % process noise covariance (nx-by-nx)
Rn = eye(ny);   % measurement noise covariance (ny-by-ny)
kest = kalman(sys_aug, Qn, Rn);

%% --- Form the LQG Tracking Controller ---
% Combine the estimator and state-feedback gain.
trksys = lqgtrack(kest, K);
% Specify the controller’s input names:
%  - First input: reference signal "r"
%  - Second input: measured plant output "y"
trksys.InputName = {'r','y'};
% (Do not reassign the OutputName property here; use the default names.)
%
% By default, lqgtrack returns a controller whose outputs are grouped as follows:
%   - The first output (in group "OutputEstimate") is the control action.
%   - The remaining outputs (in group "StateEstimate") are the observer states.
% We can retrieve these default output names as follows:
ctrlOutNames = trksys.OutputName;  % For example: {'OutputEstimate','StateEstimate1','StateEstimate2','StateEstimate3'}

%% --- Close the Loop ---
% Build the closed-loop system by connecting the plant and controller.
% The interconnection is:
%   - The controller’s output (its first output, i.e. the control action) drives the plant.
%   - The plant output "y" is fed back to the controller.
% The free external input is the reference "r".
%
% We specify external outputs as the plant’s "y" and all outputs from the controller.
clsys = connect(sys, trksys, {'r'}, [{'y'}, ctrlOutNames]);

%% --- Simulation ---
t = 0:0.01:10;           % simulation time vector (0 to 10 seconds)
r = ones(length(t), 1);   % step reference input (magnitude 1)

% Simulate the closed-loop response.
[Y, t_out, X_out] = lsim(clsys, r, t);
%
% The external outputs of clsys are:
%   Column 1: Plant output "y"
%   Column 2: Controller output channel 1 (control action)
%   Columns 3-end: Controller output channels 2,3,... (observer state estimates)
%
% In our case, if the plant output has 1 channel and the controller has 4 outputs,
% Y will have 1 + 4 = 5 columns. We assume:
%   Y(:,1)  -> y (plant output)
%   Y(:,2)  -> control action (u)
%   Y(:,3:5)-> observer states

%% --- Plot the Results ---
figure

subplot(3,1,1)
plot(t_out, Y(:,1), 'b','LineWidth',2)
xlabel('Time (s)')
ylabel('y (Plant Output)')
title('Plant Output vs. Time')
grid on

subplot(3,1,2)
plot(t_out, Y(:,2), 'r','LineWidth',2)
xlabel('Time (s)')
ylabel('u (Control Action)')
title('Control Action vs. Time')
grid on

subplot(3,1,3)
plot(t_out, Y(:,3), 'k--','LineWidth',1.5), hold on
plot(t_out, Y(:,4), 'm--','LineWidth',1.5)
plot(t_out, Y(:,5), 'g--','LineWidth',1.5)
xlabel('Time (s)')
ylabel('Observer States')
title('Observer (Kalman Filter) State Estimates vs. Time')
legend('StateEstimate1','StateEstimate2','StateEstimate3')
grid on
