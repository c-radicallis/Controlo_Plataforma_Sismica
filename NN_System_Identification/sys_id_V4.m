close all
clear all
clc

% Deep Learning and System Identification , Ljung , Andersson, Tiels ,
% Schon

data=load("TDOF_PassFD_ElCentro_s1_70_V0.mat");


%sensa: Sensibilidade dos sensores:1-Ag;2-Ai;3-As;4-Xg;5-Xgref;6-Xig;7-Xsi;
sensa=ones(8,1);  %acerto dos sinais

sname=fieldnames(data);

time=eval(['data.',sname{1},'.X.Data']);
dt=mean(diff(time));

% ddxg g is the input absolute acceleration at the base (ground);
%   and xrg={xig xsg}^T is the vector of relative displacements of
% the DOF to the ground: xig=xi-xg; xsg=xs-xg.

xg=eval(['data.',sname{1},'.Y(5).Data']);xg=double(xg)*sensa(5);
xgr=eval(['data.',sname{1},'.Y(6).Data']);xgr=double(xgr)*sensa(6);
xig=eval(['data.',sname{1},'.Y(7).Data']);xig=double(xig)*sensa(7);
xsi=eval(['data.',sname{1},'.Y(8).Data']);xsi=double(xsi)*sensa(8);
ag=eval(['data.',sname{1},'.Y(1).Data']);ag=double(ag)*sensa(1);
ai=eval(['data.',sname{1},'.Y(2).Data']);ai=double(ai)*sensa(2);
as=eval(['data.',sname{1},'.Y(3).Data']);as=double(as)*sensa(3);

% figure
% plot(time,xg)

dimv=size(xg);
if dimv(1)==1
    time=time';xg=xg'*1e-3;xgr=xgr'*1e-3;xig=xig'*1e-3;xsi=xsi'*1e-3;ag=ag';ai=ai';as=as';
end


%% Create a discrete−time neural state−space object with 
% 6 states: xig xsi ai as
states = 4
% 2 inputs: xg , ag 
Num_I=2
U = [xg  ag];
U_train = U(1: ceil( length(time)/2 ) , :  );
U_test = U(ceil( length(time)/2 ) + 1 : end, :  );
% 4 outputs:  xig xsi  ai as
Num_O=4 
Y = [ xig  xsi  ai as];
Y_train = Y(1: ceil( length(time)/2 ) , :  );
Y_test = Y(ceil( length(time)/2 ) + 1 : end, :  );
% sample time of dt seconds.

 train_data = iddata( Y_train , U_train , dt ) 
test_data =  iddata( Y_test, U_test , dt ) 
data = iddata( Y , U , dt )


%% 3.1 Example: Black-box Modeling of Robot Arm Dynamics

% To begin, prepare model regressors. Here, we pick 2
% lags in the output variable and 8 in the input variable.

% declare model variables
vars = ['xg','ag','xig','xsi','ai','as'];
% create linear regressors
R = linearRegressor(vars,{1:4,1:4});

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