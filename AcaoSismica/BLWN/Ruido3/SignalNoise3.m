%% Geracao de Ruido Aleatorio 
%%
clear all

%sinal
dt=10^-4;   %discretizacao
tfinal=20;  %tempo final

%vetor de tempos
time=(0:dt:tfinal)';

%white noise (normally dist random noise)
sig=randn(length(time),1);

%Eliminar as componentes fc e integrar/derivar
%entradas:dt,signal,cos-taper time,cross-over freq,int/deriv,show data;
[dto,yfo]=RapaFsig2(dt,sig,2,[1 50],-2,3);

%% gravar serie de ruido com amplitude unitaria
timeo=(0:dto:(length(yfo)-1)*dto)';
Ruido=[timeo yfo];
uisave('Ruido','NoiseSignal')

%% visualisar sinal
clear all

load NoiseSignal 

%entradas: sinal;janela de tempo;sobreposicao(%);esc. lin ou log;mss 
SignalSpectra(Ruido,2048,70,2)

%Sinal
%entradas: sinal;janela de tempo;sobreposicao(%);esc. lin ou log;mss ou psd;
%SignalSpectraN(Ruido,2048,70,2,2)

%% Cosine Taper
%aplicacao do cosine taper ao inicio e fim c/tempo=tc s
serf=CosTaper([time sig],3,0); %costaper 3s




