function closed_loop_sys=compute_Optimal_Controller( AA , BB , CC , DD)
sys = ss(AA,BB,CC,DD);
nx = size(AA,1);    % Number of states
nu = size(BB,2);    % Number of control inputs (should be 1)
ny = size(CC,1);    % Number of outputs

plant_aug = ss(AA, BB,[eye(nx);CC],DD);
plant_aug.InputName = {'i_sv'};   % plant input: control signal
plant_aug.OutputName = {'Qsv' , 'Fp' , 'xT' , 'x1' , 'x2','dxT' , 'dx1' , 'dx2' , 'y_xT'};  % plant output
sumblk1 = sumblk('e = x_tgt - y_xT'); % Compute the error signal: e = r - y
integrator = tf(1,[1 0]); % The integrator integrates the tracking error.
integrator.InputName = {'e'};    % error: e = r - y
integrator.OutputName = {'xi'};  % integrated error

Q = 1e3*diag([zeros(1,nx),1]);%blkdiag(eye(nx), eye(ny));
R = 1e-9*eye(nu);
K_lqi = lqi(sys, Q, R)% Design the LQI controller for the original system
K  = K_lqi(1:nx);      % state feedback gains
Ki = K_lqi(end);        % integrator gain
controller = ss([], [], [], -[K, Ki]); %   u = -[K  Ki] * [x; xi]
controller.InputName = {'Qsv' , 'Fp' , 'xT' , 'x1' , 'x2','dxT' , 'dx1' , 'dx2' , 'xi'};
controller.OutputName = {'i_sv'};
 closed_loop_sys = connect(plant_aug,  controller , integrator, sumblk1, 'x_tgt', 'y_xT')
end