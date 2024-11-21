% Deep Learning of Dynamic Systems using System Identification Toolbox™
% Dai
%%

% NeuralSS models are encapsulated by the idNeuralStateSpace
% objects. One creates a template model by specifying the
% number of states, and optionally, the number of inputs,
% number of outputs, sample time and other attributes in
% the idNeuralStateSpace constructor.

% Create a discrete−time neural state−space
% object with 3 states, 2 inputs, 4 outputs,
% and sample time of 0.1 seconds.
nss = idNeuralStateSpace(3,NumInputs=2,NumOutputs=4,Ts=0.1)

% The StateNetwork and the OutputNetwork properties of nss
% contain the deep networks representing the state-transition
% and the output functions respectively. By default, they
% employ 2 hidden layers containing 64 tanh activations
% each. The helper function createMLPNetwork facilitates easy
% creation of most commonly used networks. For example,
net = createMLPNetwork(nss,"state", LayerSizes=[4 8 4], Activations="sigmoid");
nss.StateNetwork = net;
% creates a 3-hidden-layer, sigmoid activation based network
% net and assigns it to the StateNetwork property of the model,
% which represents the function f (·) of Equation (1).

%The model can be trained using the nlssest command
nss = nlssest(z, nss, opt)



%% 2.1 Example: Black-Box Model of SI Engine Torque Dynamics
% Data for 4 input variables (Throttle position, Wastegate
% valve, Engine speed, Spark timing) and 1 state/output
% variable (Engine torque)

% To train the model, first create a discrete-time neural
% state-space object with 1 state, 4 inputs, 1 output and
% a state network f(·) with 2 layers each of size 128.
nx = 1; % number of states = number of outputs
nssModel = idNeuralStateSpace(nx,NumInputs=4);

nssModel.StateNetwork = createMLPNetwork(nssModel,"state", ...
    LayerSizes=[128 128], WeightsInitializer="glorot", ...
    BiasInitializer="zeros", Activations="tanh")


% Then segment the training data (eData object that has
% been downsampled and normalized) into multiple data
% experiments to reduce the prediction horizon and improve
% the training speed. Such segmentations can sometimes lead
% to more generalizable results.
expSize = 20;
Expts = segmentData(eData,expSize);


%Next, set up the training options of algorithm.
StateOpt = nssTrainingOptions("adam");
StateOpt.MaxEpochs = 90;
StateOpt.InputInterSample = "pchip";

%Finally, train the model using nlssest command.
nssModel = nlssest(Expts,nssModel,StateOpt)



%% 2.2 Example: Feature Reduction Using Auto-Encoders
%we assume that the true model order is unknown and choose 
% a rich set of regressors composed of the lagged I/O variables
% up to a maximum lag of 10

%In the neural
% state-space structure, we enable the use of an auto-encoder
% by setting the value of the LatentDim property to a finite
% number; this value indicates the dimension of the latent
% space.

nx = 20; % measured number of states
nu = 1; % number of inputs
nd = 7; % actual order (latent layer dim)
sys = idNeuralStateSpace(nx, NumInputs=nu, LatentDim=nd);
net1 = createMLPNetwork(sys, "state", LayerSizes=[], ...
    Activations="sigmoid", WeightsInitializer="zeros");
net2 = createMLPNetwork(sys, "encoder", LayerSizes=10, Activations="tanh");
net3 = createMLPNetwork(sys, "decoder",LayerSizes=10, Activations="tanh");
sys = setStateNetwork(sys,net1);
sys.Encoder = net2;
sys.Decoder = net3;

% Next, set up the training options.
opt = nssTrainingOptions("adam");
opt.LearnRate = 0.005;
opt.MaxEpochs = 1000;
opt.LossFcn = "MeanSquaredError";

%Finally, use the nlssest command to train the model.
sys = nlssest(data,sys,opt);
% The model sys uses 7 states. Figure 6 shows the performance
% of the model on the validation dataset.


%% 3.1 Example: Black-box Modeling of Robot Arm Dynamics

% To begin, prepare model regressors. Here, we pick 2
% lags in the output variable and 8 in the input variable.

% declare model variables
vars = ["Angular Velocity", "Torque"];
% create linear regressors
R = linearRegressor(vars,{1:2,1:8});

% Create a neural network regression function that
% uses 2 hidden layers, each containing 5 units of relU
% activation functions.
Fcn = idNeuralNetwork([5 5], "relu");

% Use the regressor set R, and the nonlinear map Fcn to
% instantiate a Nonlinear ARX model structure.
sys1 = idnlarx(vars(1),vars(2),R,Fcn);
% The weights and biases of the fully connected layers of
% Fcn constitute the unknown parameters of the model
% sys1. We want to estimate them so that the model
% response matches the training data output signal.

% Prepare training options. Use Levenberg–Marquardt
% search method and choose to minimize the simulation
% errors. This amounts to training the model in a recurrent
% setup. Also, use "zscore" as the normalization
% method for the model’s regressors and output.
opt = nlarxOptions(SearchMethod="lm",...
Focus="simulation",Display="on");
opt.NormalizationOptions.NormalizationMethod="zscore";

% Finally, train the model using the nlarx command. For
% training, split the estimation data into experiments of
% 500 samples with no overlap.
FS = 500; % data frame size
FR = FS; % frame rate
eDataSplit = segmentData(eData, FS, FR);
sys1 = nlarx(eDataSplit, sys1, opt);