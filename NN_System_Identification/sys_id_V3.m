close all
%clear all
clc

% Deep Learning and System Identification , Ljung , Andersson, Tiels ,
% Schon

addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\NN_System_Identification\dados_elcentro'
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

figure
 plot(time,xg)

dimv=size(xg);
if dimv(1)==1
    time=time';
    xg=xg'*1e-3;
    xgr=xgr'*1e-3;
    xig=xig'*1e-3;
    xsi=xsi'*1e-3;
    ag=ag';
    ai=ai';
    as=as';
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
test_data.Tstart=  time(ceil( length(time)/2 ))
figure
plot(train_data)
hold on 
plot(test_data)
legend('train','test')

test_data.Tstart=  time(ceil( length(time)/2 ))
data = iddata( Y , U , dt )


%% Frequency content

% % Assuming 'data' is your iddata object
% inputSignal = data.u; % Input signal
% outputSignal = data.y; % Output signal
% Fs = 1 / data.Ts; % Sampling frequency
% N = length(inputSignal); % Number of samples
% 
% % Compute FFT
% inputFFT = fft(inputSignal);
% outputFFT = fft(outputSignal);
% frequencies = (0:N-1)' * (Fs / N); % Frequency vector
% 
% % Only keep the first half of the spectrum (positive frequencies)
% halfIndex = floor(N/2)+1;
% frequencies = frequencies(1:halfIndex);
% inputFFT = abs(inputFFT(1:halfIndex)); % Magnitude of FFT
% outputFFT = abs(outputFFT(1:halfIndex));
% 
% % Plot the frequency content
% figure;
% subplot(2, 1, 1);
% semilogx(frequencies, inputFFT);
% title('Frequency Content of Input Signal');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% 
% subplot(2, 1, 2);
% semilogx(frequencies, outputFFT);
% title('Frequency Content of Output Signal');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');

%%  Power Spectral Density
% 
inputSignal = data.u; % Input signal
outputSignal = data.y; % Output signal
Fs = 1 / data.Ts; % Sampling frequency

% Plot PSD for input signal
figure;
subplot(2, 1, 1);
pspectrum(inputSignal, Fs, 'power');
title('Power Spectral Density of Input Signals');

% Plot PSD for output signal
subplot(2, 1, 2);
pspectrum(outputSignal, Fs, 'power');
title('Power Spectral Density of Output Signals');


%% Deep Cascade Network
net=cascadeforwardnet([6,6,6,6,6,6])
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


%%
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
