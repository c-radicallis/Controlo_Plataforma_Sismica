%Geracao de Ruido Aleatorio 

clear all

%sinal
dt=10^-3;   %discretizacao
tfinal=200;  %tempo final

%vetor de tempos
time=(0:dt:tfinal)';

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


%Vizualizacao do espetro do Sinal Filtrado
serief=[time serf output output2];
%entradas: sinal;janela de tempo;sobreposicao(%);esc. lin ou log;mss ou psd;
SignalSpectraN(serief,2048,70,1,1)
legend('Ruido','Ruido filtrado','Ruido Filtrado + CosTaper')

%gravar serie de ruido com amplitude unitaria
Ruido=[time output2];
uisave('Ruido','NoiseSignal')

%% Ajuste da amplitude
%criterio: valor RMS=1

clear all

%carregar a serie de ruido gerada para acertar a amplitude
%por forma a que o valor RMS=1

%load NoiseT200sBL50Hz       %RMS=1
%load NoiseRMS05T200sBL50Hz  %RMS=0.5
%load NoiseRMS01T200sBL50Hz   %RMS=0.1

%valor de amplitude a acertar
amp=0.317;  %3.5163;

Signal=[Ruido(:,1) amp*Ruido(:,2)];
SignalSpectra(Signal,2048,70,1)
%SignalSpectraN(Ruido2,2048,70,1,1)

%gravar
%uisave('Signal','Noise50')

%% Visualizar sinal


clear all

%carregar a serie de ruido gerada para acertar a amplitude
%por forma a que o valor RMS=1

%load NoiseT100sBL50Hz  %Band-limited 50 Hz; duracao=100s; dt=0.001s;
%load NoiseT20sBL50Hz          %Band-limited 50 Hz; duracao=20s; dt=0.001s;
%load NoiseT200sBL50Hz          %Band-limited 50 Hz; duracao=20s; dt=0.001s;
%
%
%load NoiseRMS05T200sBL50Hz      %Band-limited 50 Hz; duracao=20s; dt=0.001s; RMS=0.5
load NoiseRMS01T200sBL50Hz      %Band-limited 50 Hz; duracao=20s; dt=0.001s; RMS=0.1

SignalSpectra(Signal,2048,70,1)

