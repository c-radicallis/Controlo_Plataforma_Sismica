% Deep Learning and System Identification , Ljung , Andersson, Tiels ,
% Schon

% https://www.nonlinearbenchmark.org/#Silverbox


load SNLS80mV.mat
% V1 is input, V2 is output

fs=1e7/2^14;
t_V1V2 = [0:1/fs:(length(V1)-1)/fs];

figure(1)
plot(t_V1V2, V2,'DisplayName','output') 
hold on
plot(t_V1V2, V1,'DisplayName','input')
hold off
legend()


edat=iddata( V2(30550:38500)' , V1(30550:38500)',1/fs);

load Schroeder80mV
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

%% Linear model
mbj=bj(edat,[4 4 2 2 0]);  % Estimate Box-Jenkins polynomial model using time-domain data
[ys,Fit]=compare(vdat,mbj)
figure
resid(vdat,mbj)

%% Deep Cascade Network
net=cascadeforwardnet([6,6,6,6,6,6]);
N2=neuralnet(net);
mN2=nlarx(edat,[4 4 0],N2);
[Ys,Fit]=compare(vdat,mN2)
figure
resid(vdat,mN2)

%% LSTM Network
numHiddenUnits=8;featureDimension=1;%u(t)
layers=[sequenceInputLayer(featureDimension), lstmLayer(numHiddenUnits, 'OutputMode','sequence'), fullyConnectedLayer(1), regressionLayer];
options=trainingOptions( 'adam','MaxEpochs',1000, 'InitialLearnRate',0.01, 'Plots','training-progress');
XTrain=cell(1,1);
YTrain=cell(1,1);
XTrain{1}=edat.u';
YTrain{1}=edat.y';
net=trainNetwork(XTrain,YTrain,layers,options);


