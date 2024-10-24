% Deep Learning and System Identification , Ljung , Andersson, Tiels ,
% Schon


load elcentro.txt

%%
t = elcentro(:,1);
ddx = elcentro(:,2);
ddy = elcentro(:,3);
ts=t(2)
fs=1/ts

%%
figure(1)
plot(t, dxx,'DisplayName','ddx')
hold on
plot(t , ddy,'DisplayName','ddy')
hold off
legend()


edat=iddata( ddx(:3000)' , V1(:3000)',ts);



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