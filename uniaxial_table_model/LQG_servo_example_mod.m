clear; clc
% Design an LQG Servo Controller  - Example 
% https://www.mathworks.com/help/control/getstart/design-an-lqg-servo-controller.html

A = [0 1 0;
        0 0 1;
        1 0 0];    
nx = size(A,1);    %Number of states
B = [0.3; 1 ; -0.3];
% [0.3 1;
%         0 1;
%       -0.3 0.9]; 
nu = size(B,2) ; % Number of inputs
G = [-0.7 1.12;
        -1.17 1; 
        .14 1.5]; %0.1*ones(nx,nu);
C = [1.9 1.3 1];  
ny = size(C,1);    %Number of outputs
D = 0 ;%[0.53 -0.61];
H = [-1.2 -0.89];
%sys = ss(A,[B G],C,[D H]);
sys = ss(A,B,C,D);

obs = obsv(A,C);
r_obsv = rank(obs)
controlability = ctrb(A,B);
r_controlability = rank(controlability)

%Construct the optimal state-feedback gain using the given cost function by typing the following commands:

Q = blkdiag(eye(nx),eye(ny));
R = eye(nu);
K = lqi(ss(A,B,C,D),Q,R);

%Construct the Kalman state estimator using the given noise covariance data by typing the following commands:
Qn = eye(nu) 
Rn = eye(nu);
kest = kalman(sys,Qn,Rn)

%Connect the Kalman state estimator and the optimal state-feedback gain to form the LQG servo controller by typing the following command:
trksys = lqgtrack(kest,K)

