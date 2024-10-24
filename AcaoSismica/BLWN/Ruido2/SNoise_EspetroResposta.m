% Acoes

clear all, %close all

%----Acao de entrada-----------------------------------------------------
%escalonamento da acao
sc=1;
%
%------------10 acoes------------------------
%
%10 acoes de ruido branco gaussiano: T=50s, RMS=1, dt=1ms
% load Noise10seriesT50sBL50Hz
% In=Signal;In(:,2:end)=Signal(:,2:end)*sc;
%
%% Visualizacao de delocamentos velocidades e aceleracoes
%entradas: 1-dados; 2-correcao do delocamento; 3-gravar dados;
Integra(In,0,0)
%
%
%% Visualizacao do espectro
%entradas: 1-dados; 2-nº pontos da janela; 3-sobreposicao das janelas (%)
%          4-Escala:linear ou log;5-MSS ou PSD;
SignalSpectraN(In,2048,50,1,1)
%
%
%% Filtragem - filtrar sinais
% 
%filtragem
%entrada: 1-sinais; 2-gravar dados
%In2=FilterSignal(In,0);
%
In2=FilterSignal2(In,0);
%
%Vizualizacao do espetro do Sinal Filtrado
SignalSpectraN(In2,4096,70,2)

%% Comparacao dos espetros
SignalSpectraN([In In2(:,2:end)],4096,70,2)


%% Comparacado dos sinais
%entradas: Estrutura de dados com os sinais:
%
se=1; %serie a visualizar
%
Infe(1).data=[In(:,1),In(:,se+1)];
Infe(2).data=[In2(:,1),In2(:,se+1)];
%
CompSinal(Infe)

%% Espetros de Resposta

%Visualizacao do espetro
%entradas: 1-dados; 2-gravar dados;
global m k c
ResponseSpectraBuilt(In,1);

%% Visualizar Espetros de Resposta
%%
clear all

%10 acoes do tipo 1
load RespSpectra_Noise10seriesT50sBL50Hz

%estradas:espetro de resposta;comparacao c/EC8 spectra; acao tipo;
ResponseSpectraRead(data,0,1)

