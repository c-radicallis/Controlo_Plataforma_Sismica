%Geracao de Ruido Aleatorio 

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%definicao do sinal

%nº de series
nseries=10;

%valor de amplitude a acertar
amp=3.3;

%sinal
dt=10^-3;   %discretizacao
tfinal=50;  %tempo final

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%vetor de tempos
time=(0:dt:tfinal)';

Signal=zeros(length(time),nseries);
Signal(:,1)=time;

for i=1:nseries
%white noise (normally dist random noise)
sig=randn(length(time),1);

%aplicacao do cosine taper ao inicio e fim c/tempo=tc s
serf=CosTaper([time sig],3,0); %costaper 3s

%figure,plot(time,sigc)

%Espetro do Sinal
%entradas: sinal;janela de tempo;sobreposicao(%);esc. lin ou log;mss ou psd;
%SignalSpectraN([time serf],2048,70,1,2)


%sinal filtrado a 50 Hz
input=[time serf];
sim('LP50hz')

%aplicacao do cosine taper ao inicio e fim c/tempo=tc s
output2=CosTaper([time output],0.5,0); %costaper 0.5s


if nseries==1
SignalSpectra([time amp*output2],2048,70,1)
end

Signal(:,i+1)=amp*output2;

end
%gravar serie de ruido com amplitude unitaria
uisave('Signal','NoiseSignal')


%% Visualizar sinal

clear all

%carregar a serie de ruido gerada
ser=10;

load Noise10seriesT50sBL50Hz  %Band-limited 50 Hz; duracao=100s; dt=0.001s;

SignalSpectra([Signal(:,1) Signal(:,ser+1)],2048,70,1)

