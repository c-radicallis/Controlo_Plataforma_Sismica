clear; clc
% Design an LQG Servo Controller  - Example 
% https://www.mathworks.com/help/control/getstart/design-an-lqg-servo-controller.html

% Define the original system
A = [0 1 0;
     0 0 1;
     1 0 0];    
nx = size(A,1);    % Number of states
B = [0.3; 1; -0.3];
nu = size(B,2);    % Number of control inputs (should be 1)
C = [1.9 1.3 1];  
ny = size(C,1);    % Number of outputs
D = 0;

% Define process noise matrices
G = eye(nx);         % Process noise matrix (assuming full-state noise)
H = zeros(ny, nx);   % No direct noise feedthrough

% Create the augmented system
sys_aug = ss(A, [B G], C, [D H]);

% Assign input groups:
% - Channel 1 is control,
% - Channels 2 to (nx+1) are noise.
sys_aug.InputGroup.control = 1;
sys_aug.InputGroup.noise   = 2:(nx+1);
% Explicitly set the KnownInput group to only the control input
sys_aug.InputGroup.KnownInput = 1;

% Design the LQI controller for the original system
Q = blkdiag(eye(nx), eye(ny));
R = eye(nu);
K = lqi(ss(A, B, C, D), Q, R)

% Define noise covariance data
% Here Qn should be for process noise (nx-by-nx) and Rn for measurement noise (ny-by-ny)
Qn = eye(nx);
Rn = eye(ny);

% Construct the Kalman estimator using the augmented system
kest = kalman(sys_aug, Qn, Rn)

% Now, lqgtrack expects that the number of rows in K (control actions) matches 
% the number of known input channels in kest (which is now 1).
trksys = lqgtrack(kest, K)
