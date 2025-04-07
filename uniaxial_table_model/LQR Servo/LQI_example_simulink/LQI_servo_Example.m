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
plant_aug = ss(A, B,[eye(2);C], D);
plant_aug.InputName = {'u'};
plant_aug.OutputName = {'x1','x2','y'};

%% Define the Summing Block for Error Calculation
% Compute the error signal: e = r - y
sumblk1 = sumblk('e = r - y');

%% Define the Integrator
% The integrator integrates the tracking error.
integrator = tf(1,[1 0]);
integrator.InputName = {'e'};    % error: e = r - y
integrator.OutputName = {'xi'};  % integrated error

%% Compute the LQI Gains
% Choose weighting matrices for states (including integrator) and input
Q = diag([0, 0, 1]);  % weights for [x1; x2; xi]
R = 1e-1;                % weight for input

% Compute the LQI gain: [K  Ki] where u = -K*x - Ki*xi
K_lqi = lqi( plant, Q, R)
K  = K_lqi(1:2);      % state feedback gains
Ki = K_lqi(3);        % integrator gain

% Define the Controller Block
% The controller is a static gain that uses the plant states and the integrator state:
%   u = -[K  Ki] * [x; xi]
controller = ss([], [], [], -[K, Ki]);
controller.InputName = {'x1','x2','xi'};
controller.OutputName = {'u'};

%% Interconnect the Blocks Using connect
% The closed-loop interconnection consists of:
% - sumblk1:   r - y => e
% - integrator: e => xi
% - controller:  -K*x - Ki*xi => u
% - plant_aug: u => x1 , x2 , y

% Here, 'r' is the external reference input and 'y' is the measured output.
clsys = connect(plant_aug,  controller , integrator, sumblk1, 'r', 'y')

%% Simulate the Closed-Loop Response
t = 0:0.005:5;          % time vector from 0 to 10 seconds
r = ones(length(t), 1); % step reference

% Simulate the closed-loop response using lsim:
[y_out, t_out, x] = lsim(clsys, r, t);

% Plot the response:
figure(1)
hold on
plot(t,r,'-.')
plot(t_out,x(:,1:2))
plot(t_out,-K_lqi*x','--') % u = -K_lqi*x'
xlabel('Time (s)')
ylabel('Output y')
title('Closed-Loop Response with LQG Tracking Controller')
grid on
legend

%% % %checking if connect is correct
% K_sys = connect(controller,plant_aug,'xi','y')
% 
% K_sys.A == (A-B*K)
% K_sys.B == (-B*Ki)