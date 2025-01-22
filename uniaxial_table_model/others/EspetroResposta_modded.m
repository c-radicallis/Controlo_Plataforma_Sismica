%function EspetroResposta_function( nome_acao )
% o scrip 'EspetroResposta18102024.m' irá permitir-te carregar o sismo com a função 'load sismo.txt', e no final tens a funcao 'ResponseSpectra()'
% que te permite visualizar o espetro de resposta e gravar. O file gravado pode depois ser visualizado conforme está no final do script. Está ainda um passo 
% intermedio com a funcao 'Integra()' que te permite visualizar os deslocamentos, velocidades e aceleracoes (que podes ou não correr).

%Script para construir/visualizar o espetro de resposta

%clear all, %close all

% ------------Inputs---------------------------------------------
% carregar file do sismo: load sismo.txt -> criada matriz 'sismo' com dados 
% dados: 1ªcoluna-vetor de tempos; 2ª a nª coluna-aceleracoes
% fazer In=sismo para facilitar na utilizacao dos sripts

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
%if nome_acao == 'elcentro'
%Sismo Elcentro
%cl=3; %2-FN;3-FP
load elcentro.txt
In=elcentro;
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
%cl=2; %2-FN;3-FP
%load sylmar.txt
%In=sylmar;
% In(:,1)=sylmar(:,1);In(:,2)=sylmar(:,cl)*sc;
%
%% Visualizacao de delocamentos velocidades e aceleracoes
%entradas: 1-dados; 2-correcao do delocamento; 3-gravar dados;
%dados: In-matriz de dados: 1ªcoluna-vetor de tempos; 2ª a nª coluna-aceleracoe


%Integra(In,0,0)
%
%% Espetros de Resposta
%entradas: 1-dados; 2-gravar dados: 0-nao, 1-sim;
%dados: In-matriz de dados: 1ªcoluna-vetor de tempos; 2ª a nª coluna-aceleracoe

% global m k c
% ResponseSpectraBuilt(In,0);

%% Visualizar Espetros de Resposta Gravados

%clear all

% ------------Inputs---------------------------------------------
%visualizar espetro medio: 0-nao, 1-sim
spm=1;

% 'file.mat' gravado em 'ResponseSpectraBuilt(In,1)':
load ElCentro_RespSpectra

% -----------------------------------------------------------------

%Espetro de Resposta medio em Deslocamento,Velocidade e Aceleracao

if spm==1
figure('Position',[50 50 1200 400]),
subplot(121),semilogx(data.T,data.mdr),ylabel('D (m)'),xlabel('T(s)'),title('medium response spectra')
subplot(122),semilogx(data.T,data.m),ylabel('A (m/s^2)'),xlabel('T(s)'),title('medium response spectra')
end

figure('Position',[50 50 1200 400]),
subplot(121),semilogx(data.T,data.dr),ylabel('D (m)'),xlabel('T(s)'),title('response spectra')
subplot(122),semilogx(data.T,data.a),ylabel('A (m/s^2)'),xlabel('T(s)'),title('response spectra')

%end