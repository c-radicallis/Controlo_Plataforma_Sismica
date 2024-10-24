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
% load EC8_A1_f.mat
% In=data;In(:,2:end)=data(:,2:end)*sc;
%
%10 acoes EC8 DNA Zona 2.1 Solo D Categ. II
% load EC8_A2.txt
% In=EC8_A2;In(:,2:end)=EC8_A2(:,2:end)*sc;
%
%
%10 acoes de ruido branco gaussiano: T=50s, RMS=1, dt=1ms
% load Noise10seriesT50sBL50Hz
% In=Signal;In(:,2:end)=Signal(:,2:end)*sc;
%
%-----------1 acao--------------------------------
%
%acao EC8 DNA Zona 1.1 Solo D Categ. II
% load EC8DNA11DIIacel.txt  %sinal sismo - valores em m/s2
% In(:,1)=EC8DNA11DIIacel(:,1);In(:,2)=EC8DNA11DIIacel(:,2)*sc;  
%
%acao EC8 DNA Zona 2.1 Solo D Categ. II
% load EC8DNA21DIIacel.txt  %sinal sismo - valores em m/s2
% In(:,1)=EC8DNA21DIIacel(:,1);In(:,2)=EC8DNA21DIIacel(:,2)*sc;
%
%acao EC8
% load SinalAcel.txt  %sinal sismo - valores em g
% In(:,1)=SinalAcel(:,1);In(:,2)=SinalAcel(:,2)*9.86*sc;
%
%Sismo Elcentro
% cl=3; %2-FN;3-FP
% load elcentro.txt
% In=elcentro;
%In(:,1)=elcentro(:,1);In(:,2)=elcentro(:,cl)*sc;
%
%Sismo Erzikan
% cl=3; %2-NS;3-EW
% load erzikan.txt
% In=erzikan;
% In(:,1)=erzikan(:,1);In(:,2)=erzikan(:,cl)*sc;
%
%Sismo Jiji
% cl=3; %2-NS;3-EW
% load jiji.txt
% In=jiji;
% In(:,1)=jiji(:,1);In(:,2)=jiji(:,cl)*sc;
%
%Sismo Kobe
% cl=3; %2-NS(FN);3-EW(FP)
% load kobe.txt
% In=kobe;
% In(:,1)=kobe(:,1);In(:,2)=kobe(:,cl)*sc;
%
%Sismo Newhall
% cl=3; %2-FN;3-FP
% load newhall.txt
% In=newhall;
% In(:,1)=newhall(:,1);In(:,2)=newhall(:,cl)*sc;
%
%Sismo Rinaldi
% cl=3; %2-FN;3-FP
% load rinaldi.txt
% In=rinaldi;
% In(:,1)=rinaldi(:,1);In(:,2)=rinaldi(:,cl)*sc;
%
%Sismo Sylmar
cl=2; %2-FN;3-FP
load sylmar.txt
In=sylmar;
% In(:,1)=sylmar(:,1);In(:,2)=sylmar(:,cl)*sc;
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

clear all

%acao tipo 1 ou 2
sa=1;

if sa==1
%10 acoes do tipo 1
load ResponseSpectraEC8_A1_f
elseif sa==2
%10 acoes do tipo 2
load ResponseSpectraEC8_A2
end

%estradas:espetro de resposta;comparacao c/EC8 spectra; acao tipo;
ResponseSpectraRead(data,1,sa)

%Espetro de Resposta em Deslocamento,Velocidade e Aceleracao
vm=data.m.*data.T'/2/pi;
dm=vm.*data.T'/2/pi;

figure('Position',[50 50 1200 400]),
subplot(131),semilogx(data.T,dm),ylabel('D (m)'),xlabel('T(s)')
subplot(132),semilogx(data.T,vm),ylabel('V (m/s)'),xlabel('T(s)')
subplot(133),semilogx(data.T,data.m),ylabel('A (m/s^2)'),xlabel('T(s)')



%% Espetros de Resposta de Sismos Comuns
%% Sismo ElCentro
clear all
data1=load('RespSpectra_ElCentroFN');
data2=load('RespSpectra_ElCentroFP');
%
%entradas: Espetros de Resposta; Intervalo de periodos
ResponseSpectraRead1([data1 data2],[0 5])
legend('FN','FP')

%% Sismo Erzikan
clear all
data1=load('RespSpectra_ErzikanNS');
data2=load('RespSpectra_ErzikanEW');
%
%entradas: Espetros de Resposta; Intervalo de periodos
ResponseSpectraRead1([data1 data2],[0 5])
legend('FN','FP')

%% Sismo Jiji
clear all
data1=load('RespSpectra_JijiNS');
data2=load('RespSpectra_JijiEW');
%
%entradas: Espetros de Resposta; Intervalo de periodos
ResponseSpectraRead1([data1 data2],[0 5])
legend('FN','FP')

%% Sismo Kobe
clear all
data1=load('RespSpectra_KobeNS');
data2=load('RespSpectra_KobeEW');
%
%entradas: Espetros de Resposta; Intervalo de periodos
ResponseSpectraRead1([data1 data2],[0 5])
legend('FN','FP')

%% Sismo Newhall
clear all
data1=load('RespSpectra_NewhallFN');
data2=load('RespSpectra_NewhallFP');
%
%entradas: Espetros de Resposta; Intervalo de periodos
ResponseSpectraRead1([data1 data2],[0 5])
legend('FN','FP')

%% Sismo Rinaldi
clear all
data1=load('RespSpectra_RinaldiFN');
data2=load('RespSpectra_RinaldiFP');
%
%entradas: Espetros de Resposta; Intervalo de periodos
ResponseSpectraRead1([data1 data2],[0 5])
legend('FN','FP')

%% Sismo Sylmar
clear all
data1=load('RespSpectra_SylmarFN');
data2=load('RespSpectra_SylmarFP');

%entradas: Espetros de Resposta; Intervalo de periodos
ResponseSpectraRead1([data1 data2],[0 5])
legend('FN','FP')
