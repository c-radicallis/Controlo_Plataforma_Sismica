close all
clear all
clc

% Deep Learning and System Identification , Ljung , Andersson, Tiels ,
% Schon

x_i=0
v_i=10
a=-9.81

t_f=-2*v_i/a
time = 0:0.001:t_f;

v=v_i+a*time;
x=x_i+v_i*time+0.5*a*time.^2;

x_noisy = awgn(x,30,"measured");
a_noisy = awgn(a*ones(size(time)),40,"measured");

figure
plot(time,x)
hold on
plot(time,x_noisy)
hold off

figure
plot(time,a)
hold on
plot(time,a_noisy)

x=x';
a=a*ones(size(time))';
x_noisy=x_noisy';
a_noisy=a_noisy';

dt=mean(diff(time))

data =iddata([x a],[x_noisy a_noisy],dt)

%% Create a discrete−time neural state−space object with 
% %states: x a
% states = 2

% inputs: x_noisy  a_noisy
Num_I=2
U = [x_noisy  a_noisy];
U_train = U(1: ceil( length(time)/2 ) , :  );
U_test = U(ceil( length(time)/2 ) + 1 : end, :  );
% outputs:  x a
Num_O=2 
Y = [ x  a ];
Y_train = Y(1: ceil( length(time)/2 ) , :  );
Y_test = Y(ceil( length(time)/2 ) + 1 : end, :  );
% sample time of dt seconds.

 train_data = iddata( Y_train , U_train , dt ) 
test_data =  iddata( Y_test, U_test , dt ) 
data = iddata( Y , U , dt )


%% Deep Cascade Network
net=cascadeforwardnet([6,6,6,6,6,6])
view(net)
N2=neuralnet(net)

%%
N_A = 10*ones(Num_O)
N_B = 10*ones(Num_O,Num_I)
N_K = zeros(Num_O,Num_I)

% opt = nlarxOptions
% opt.SearchOptions.MaxIterations = 10
mN2=nlarx(data,[N_A N_B N_K ],N2)

%% Compare
compare(test_data,mN2)

figure(4)
resid(data,mN2)


%% Validation
% array to timetable


idplot(test_data)

%% simOpt = simOptions('InitialCondition',[0 0 0 0]);
yn = sim( mN2 , test_data );

%%
figure
plot( time(ceil( length(time)/2 ):end) , Y_test(:,1))
hold on

%%
plot(yn(:,1))
xlabel("Time"); ylabel("State");
legend("Original","Estimated");




%% Compare
% 
% vdat = iddata(Y_test,U_test,dt)
% [Ys,Fit]=compare(vdat,nss)

%% 2.2 Example: Feature Reduction Using Auto-Encoders
% %we assume that the true model order is unknown and choose
% % a rich set of regressors composed of the lagged I/O variables up to a
% % maximum lag of 10
% 
% %In the neural
% % state-space structure, we enable the use of an auto-encoder by setting
% % the value of the LatentDim property to a finite number; this value
% % indicates the dimension of the latent space.
% 
% nx = 20; % measured number of states
% nu = 1; % number of inputs
% nd = 7; % actual order (latent layer dim)
% sys = idNeuralStateSpace(nx, NumInputs=nu, LatentDim=nd);
% net1 = createMLPNetwork(sys, "state", LayerSizes=[], ...
%     Activations="sigmoid", WeightsInitializer="zeros");
% net2 = createMLPNetwork(sys, "encoder", LayerSizes=10, Activations="tanh");
% net3 = createMLPNetwork(sys, "decoder",LayerSizes=10, Activations="tanh");
% sys = setStateNetwork(sys,net1);
% sys.Encoder = net2;
% sys.Decoder = net3;
% 
% % Next, set up the training options.
% opt = nssTrainingOptions("adam");
% opt.LearnRate = 0.005;
% opt.MaxEpochs = 1000;
% opt.LossFcn = "MeanSquaredError";
% 
% %Finally, use the nlssest command to train the model.
% sys = nlssest(data,sys,opt);
% % The model sys uses 7 states. Figure 6 shows the performance of the model
% % on the validation dataset.
