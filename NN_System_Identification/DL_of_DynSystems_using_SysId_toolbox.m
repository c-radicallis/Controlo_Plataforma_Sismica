% Deep Learning of Dynamic Systems using System Identification Toolbox™

%%
% Create a discrete−time neural state−space
% object with 3 states, 2 inputs, 4 outputs,
% and sample time of 0.1 seconds.
nss = idNeuralStateSpace(3,NumInputs=2,NumOutputs=4,Ts=0.1)

net = createMLPNetwork(nss,"state", LayerSizes=[4 8 4], Activations="sigmoid");
nss.StateNetwork = net;

nss = nlssest(z, nss, opt)

%%

nx = 1; % number of states = number of outputs
nssModel = idNeuralStateSpace(nx,NumInputs=4);
nssModel.StateNetwork = createMLPNetwork(nssModel,"state", LayerSizes=[128 128], WeightsInitializer="glorot", BiasInitializer="zeros", Activations="tanh")

expSize = 20;
Expts = segmentData(eData,expSize);
Next, set up the training options of algorithm.
StateOpt = nssTrainingOptions("adam");
StateOpt.MaxEpochs = 90;
StateOpt.InputInterSample = "pchip";

nssModel = nlssest(Expts,nssModel,StateOpt)