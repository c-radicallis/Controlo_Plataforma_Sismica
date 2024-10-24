% Acoes

clear all, %close all

%----Acao de entrada-----------------------------------------------------
%escalonamento da acao
sc=1;
%
%------------10 acoes------------------------
%
%10 acoes EC8 DNA Zona 1.1 Solo D Categ. II
% load EC8_A1.txt
% In=EC8_A1;In(:,2:end)=EC8_A1(:,2:end)*sc;
%
% load EC8_A1_f.mat
% In=data;In(:,2:end)=data(:,2:end)*sc;
%
%
%10 acoes EC8 DNA Zona 2.1 Solo D Categ. II
load EC8_A2.txt
In=EC8_A2;In(:,2:end)=EC8_A2(:,2:end)*sc;

%% serie a visualizar
se=1;
sinal=[In(:,1),In(:,se+1)];

clear In
In=sinal;

%% Visualizacao de delocamentos velocidades e aceleracoes
%entradas: 1-dados; 2-correcao do delocamento; 3-gravar dados;
Integra(In,0,1)
%
%
%% Visualizacao do espectro
%entradas: 1-dados; 2-nº pontos da janela; 3-sobreposicao das janelas (%)
%          4-Escala:linear ou log;5-MSS ou PSD;
SignalSpectraN(In,2048,50,2,2)
%
%% Preparar Acoes para dSPACE

clear all

%load SAtp1se4
%load SAtp1fse1
load SAtp2se1

tin=20; %offset no tempo

Ruido=[(data.time+tin)';data.x'];

uisave('Ruido','NoiseSignal')

