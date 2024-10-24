%% Geracao de Ruido Aleatorio 4
%%
clear all

%sinal
dt=10^-3;   %discretizacao
tfinal=100;  %tempo final
ct=2;       %cosine-taper time
fc=[0.3 1000]; %intervalo de frequencias
di=-2;       %integrar/derivar o sinal: -2,-1,0,1,2


%Parametros
L = tfinal/dt+1;                     % Length of signal
NFFT = 2^nextpow2(L)*2+1;                % Next power of 2 from length of y

%Espetro do Sinal
%sinal de Amplitude Cte
f = 1/dt*linspace(0,1,NFFT); %vetor de frequencias
Amp=ones(NFFT,1);
Pha=pi*randn(NFFT,1);
yf=Amp.*exp(complex(zeros(NFFT,1),Pha));  %espetro

%sinal no tempo
yt=ifft(yf,'symmetric')*L;


df=f(2)-f(1);           %discretizacao em frequencia
nfc1=fix(fc(1)/df)+1;   %pt da freq 1
nfc2=ceil(fc(2)/df)+1;  %pt da freq 2

%eliminacao dos termos fora da banda de interesse fc
yf2=[zeros(nfc1-1,1);yf(nfc1:nfc2);zeros(NFFT-nfc2,1)];

%sinal no tempo
yt2=ifft(yf2,'symmetric')*L;


%aplicacao do integral/derivada
wf=2*pi*f(2:end); %freq em radianos

%TF
yf03=(complex(zeros(size(f(2:end)')),wf')).^di.*yf2(2:end);
yf3=[0;yf03];


%serie temporal
yt03 = ifft(yf3,'symmetric')*L;


%tempo
time=(0:dt:(L-1)*dt)';

%aplicacao do cosine taper ao inicio e fim c/tempo=tc s
yt3=CosTaper([time yt03(1:length(time))],ct,0); 


figure('Position',[90 350 600 300])
subplot(221),semilogy(f,2*abs(yf),f(nfc1:nfc2),2*abs(yf(nfc1:nfc2))),xlim([0 2*fc(2)])
title('random'),legend('original','filter') 
subplot(222),semilogy(f,2*abs(yf3)),title('random+int/der'),xlim([0 2*fc(2)])
subplot(223),plot(f,angle(yf),f(nfc1:nfc2),angle(yf(nfc1:nfc2))),xlim([0 2*fc(2)])
subplot(224),plot(f,angle(yf3)),xlabel('f(Hz)'),xlim([0 2*fc(2)]),xlabel('f(Hz)'),

figure('Position',[700 50 600 600]),
subplot(311),plot(time,yt(1:length(time))), title('random')
subplot(312),plot(time,yt2(1:length(time))), title('random+filter')
subplot(313),plot(time,yt3), title('random+filter+int/der')
xlabel('t(s)')

%% gravar serie de ruido com amplitude unitaria

tin=20;

Ruido=[(time+tin)';yt3'];

uisave('Ruido','NoiseSignal')

%% visualisar sinal
clear all

%load NoiseSignal 
%load NoiseSignalFpa03hzT100s 
%load NoiseSignalFpa1hzT100s 
load NoiseSignalFpa075hzT100s 


%entradas: sinal;janela de tempo;sobreposicao(%);esc. lin ou log;mss 
SignalSpectra(Ruido',2048,50,2)

%Sinal
%entradas: sinal;janela de tempo;sobreposicao(%);esc. lin ou log;mss ou psd;
%SignalSpectraN(Ruido,2048,70,2,2)

%%
%Eliminar as componentes fc e integrar/derivar
%entradas:dt,signal,cos-taper time,cross-over freq,int/deriv,show data;
clear
load NoiseSignalFpa03hzT100s

dt=mean(diff(Ruido(1,:)));
[dto,yfo]=RapaFsig2(dt,Ruido(2,:)',0.1,[0.75 500],0,1);

tf=0:dto:Ruido(1,end)-Ruido(1,1);
serf=CosTaper([tf' yfo(1:length(tf))],1,0); %costaper 3s


tin=20;
tf2=tf'+tin;

figure,plot(Ruido(1,:),Ruido(2,:),tf2,serf)

%% gravar o ruido filtrado
clear Ruido
Ruido=[tf2';serf'];
uisave('Ruido','NoiseSignal')



%% Espetro do Sinal
%entradas: sinal;janela de tempo;sobreposicao(%);esc. lin ou log;mss 
SignalSpectra([time+tin yt3],2048,70,2)

%comparacao com o original (em aceleracao)
SignalSpectra([time(1:end-2) diff(yt3,2)],2048,70,2)

%no tempo
figure,plot(time,yt2(1:length(time)),time(1:length(time)-2),diff(yt3,2)/dt^2)
xlabel('t(s)')



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




