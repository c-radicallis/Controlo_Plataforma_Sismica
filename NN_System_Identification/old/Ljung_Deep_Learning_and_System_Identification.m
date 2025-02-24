clear all

% Deep Learning and System Identification 
% Ljung , Andersson, Tiels , Schon

% https://www.nonlinearbenchmark.org/#Silverbox

% README
% 
% Additional data of the Silverbox system, not described in the technical
% note (http://www.it.uu.se/research/publications/reports/2013-006/2013-006-nc.pdf)
% can be found in Schroeder80mV.mat or Schroeder80mV.csv. This dataset
% contains, amongst others, a Schroeder phase multisine measurement. This
% one period of the schroeder phase multisine can be extracted with the
% following (Matlab code). The measurement settings are the same as for the
% data described in the technical note.

load('Schroeder80mV.mat')

V1=V1-mean(V1); % Remove offset errors on the input measurements (these are visible in the zero sections of the input)
                % The input is designed to have zero mean
V2=V2-mean(V2); % Approximately remove the offset errors on the output measurements. 
                % This is an approximation because the silverbox can create itself also a small DC level 
 
uSchroeder=V1(10585:10585+1023);  % select the Schroeder section of the experiment
ySchroeder=V2(10585:10585+1023);
% One period is 1024 points. Only the odd frequencies bins (f0,3f0,5f0,...) 
% are excited. f0 = fs/N, N=1024.


load SilverBoxFiles/SNLS80mV.mat
% V1 is input, V2 is output

fs=1e7/2^14;
t_V1V2 = [0:1/fs:(length(V1)-1)/fs];

% figure(1)
% plot(t_V1V2, V2,'DisplayName','output') 
% hold on
% plot(t_V1V2, V1,'DisplayName','input')
% hold off
% legend()


edat=iddata( V2(30550:38500)' , V1(30550:38500)',1/fs);

load SilverBoxFiles/Schroeder80mV
t_Schroeder = [0:1/fs:(length(V1)-1)/fs];

figure(2)
plot(t_Schroeder,V2,'DisplayName','output')
hold on
plot(t_Schroeder,V1,'DisplayName','input')
hold off
legend()

ySchroeder=V2(10585:10585+1023);
uSchroeder=V1(10585:10585+1023);
vdat=iddata(ySchroeder' , uSchroeder' ,1/fs);

% %% Linear model
% 
% % Estimate Box-Jenkins polynomial model using time-domain data
% mbj=bj(edat,[4 4 2 2 0]);  
% [ys,Fit]=compare(vdat,mbj)
% figure(3)
% resid(vdat,mbj)

%% Deep Cascade Network
net=cascadeforwardnet([6,6,6,6,6,6]);
N2=neuralnet(net);
mmN2=nlarx(edat,[4 4 0],N2);

%%
compare(vdat,mmN2)
figure(4)
resid(vdat,mmN2)

%% LSTM Network
% numHiddenUnits=8;
% featureDimension=1;%u(t)
% layers=[sequenceInputLayer(featureDimension), lstmLayer(numHiddenUnits, 'OutputMode','sequence'), fullyConnectedLayer(1), regressionLayer];
% options=trainingOptions( 'adam','MaxEpochs',200, 'InitialLearnRate',0.01, 'Plots','training-progress');
% XTrain=cell(1,1);
% YTrain=cell(1,1);
% XTrain{1}=edat.u';
% YTrain{1}=edat.y';
% net=trainNetwork(XTrain,YTrain,layers,options);

%%
lstm_predict = predict(net,vdat.InputData');

Ts = vdat.Ts; % Sampling time
N = length(vdat.y); % Number of data points (assuming single-output data)
timeVector = (0:N-1)' * Ts; % Create the time vector

plot(timeVector,lstm_predict)
hold on
plot(timeVector,vdat.OutputData)
legend('LSTM predicted','Validation Data')