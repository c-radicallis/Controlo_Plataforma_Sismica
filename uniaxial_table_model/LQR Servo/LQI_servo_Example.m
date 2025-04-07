%% LQI Closed-Loop Design Using connect and lqi
clc; clear; close all;

%% Define the Plant
% Example plant: second-order system
A = [0 1; -2 -3];
B = [0; 1];
C = [1 0];
D = 0;
plant = ss(A, B, C, D);

% For interconnection, augment the plant to output its states.
% This allows the controller to access the full state vector.
plant_aug = ss(A, B, [C; eye(2)], D);
plant_aug.InputName = {'u'};
plant_aug.OutputName = {'x1','x2','y'};
plant_aug

%% Define the Integrator
% The integrator integrates the tracking error.
integrator = tf(1,[1 0]);
integrator.InputName = {'e'};    % error: e = r - y
integrator.OutputName = {'xi'};  % integrated error
integrator

%% Compute the LQI Gains
% Choose weighting matrices for states (including integrator) and input
Q = diag([0, 0, 1]);  % weights for [x1; x2; xi]
R = 1;                % weight for input

% Compute the LQI gain: [K  Ki] where u = -K*x - Ki*xi
K_lqi = lqi( plant, Q, R)
K  = K_lqi(1:2);      % state feedback gains
Ki = K_lqi(3);        % integrator gain

%% Define the Controller Block
% The controller is a static gain that uses the plant states and the integrator state:
%   u = -[K  Ki] * [x; xi]
controller = ss([], [], [], -[K, Ki]);
controller.InputName = {'x1','x2','xi'};
controller.OutputName = {'u'};

%% Define the Summing Block for Error Calculation
% Compute the error signal: e = r - y
sumblk1 = sumblk('e = r - y');

%% Interconnect the Blocks Using connect
% The closed-loop interconnection consists of:
% - plant_aug: Input 'u'; outputs 'y', 'x1', and 'x2'
% - integrator: Input 'e'; output 'xi'
% - controller: Inputs 'x1', 'x2', 'xi'; output 'u'
% - sumblk1: Forms error 'e' from reference 'r' and plant output 'y'
%
% Here, 'r' is the external reference input and 'y' is the measured output.
clsys = connect(plant_aug, integrator, controller, sumblk1, 'r', 'y');

%% Simulate the Closed-Loop Response
figure;
step(clsys);
title('Closed-Loop Step Response with LQI Controller');
grid on;

%% Display the Closed-Loop System
disp('Closed-loop system:');
clsys

%%
K_sys = connect(controller,plant_aug,{'x1','x2','xi'},'y')
