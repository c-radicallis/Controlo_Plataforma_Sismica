close all;clear;clc

A = [0 1 0;
     0 0 1;
     1 0 0];    
B = [0.3; 1; -0.3];
C = [0 1 0];  
D = 0;
sys = ss(A, B, C, D);

nx = size(A,1);  
ny = size(C,1);

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
% Define process noise matrices
G = eye(nx);         % Process noise matrix (assumed affecting all states)
H = zeros(ny, nx);   % No direct noise feedthrough

Q = blkdiag(eye(nx), eye(ny));
R = eye(size(B,2));  % Note: size(B,2) is the number of control inputs (should be 1)
K_lqr = lqr(A, [B G],Q,R);

