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
%detrend()

% figure
%  plot(time,xg)

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

%%
% MATLAB script to splice the beginning and end of a signal where values are almost constant

% Example input signal (can be replaced with real data)
signal = xg ;
n = length(signal);
t = 1:n; % Time or index vector

% Parameters to detect constant regions
window_size = 10; % Size of the moving window
threshold = 0.0001; % Threshold for detecting constant values (tolerance for variation)

% Calculate the local standard deviation of the signal
local_std = movstd(signal, window_size);

% Identify constant regions (where local standard deviation is below the threshold)
constant_regions = local_std < threshold;

% Find the indices for the beginning and end constant regions
start_index = find(~constant_regions, 1, 'first'); % First non-constant index
end_index = find(~constant_regions, 1, 'last'); % Last non-constant index

% Splice the signal to remove the constant regions at the beginning and end
spliced_signal = signal(start_index:end_index);

% Plot the results
figure;
subplot(3, 1, 1);
plot(t, signal);
title('Original Signal');
xlabel('Time'); ylabel('Amplitude');

yline(signal(1), '--r', 'Start');
yline(signal(end), '--r', 'End');

subplot(3, 1, 2);
plot(t, local_std);
title('Local Standard Deviation');
xlabel('Time'); ylabel('Std Dev');
yline(threshold, '--r', 'Threshold');

subplot(3, 1, 3);
new_t = t(start_index:end_index);
plot(new_t, spliced_signal);
title('Spliced Signal');
xlabel('Time'); ylabel('Amplitude');

% Display indices of the spliced regions
fprintf('The signal was spliced from index %d to %d.\n', start_index, end_index);



%% Create a discrete−time neural state−space object with 
% 6 states: xig xsi ai as

% 2 inputs: xg , ag 
Num_I=2
U = [xg(start_index:end_index)  ag(start_index:end_index)];
% 4 outputs:  xig xsi  ai as
Num_O=1 
Y = [ai(start_index:end_index)];
% sample time of dt seconds.
data = iddata( Y , U , dt )



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
figure
compare(data,mN2)


%%
figure
resid(data,mN2)


%% simOpt = simOptions('InitialCondition',[0 0 0 0]);
yn = sim( mN2 , data );

figure
plot( time(ceil( length(time)/2 ):end) , Y_test(:,1))
hold on
plot(yn(:,1))
xlabel("Time"); ylabel("State");
legend("Original","Estimated");



