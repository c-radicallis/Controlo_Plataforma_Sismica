clear;
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\AcaoSismica\Sismos'

%----Acao de entrada-----------------------------------------------------
%escalonamento da acao
sc=1;

%------------10 acoes------------------------

%% 10 acoes EC8 DNA Zona 1.1 Solo D Categ. II
load EC8_A1.txt
In=EC8_A1;In(:,2:end)=EC8_A1(:,2:end)*sc;
load EC8_A1_f.mat
In=data;In(:,2:end)=data(:,2:end)*sc;

%% 10 acoes EC8 DNA Zona 2.1 Solo D Categ. II
load EC8_A2.txt
In=EC8_A2;In(:,2:end)=EC8_A2(:,2:end)*sc;


%% 10 acoes de ruido branco gaussiano: T=50s, RMS=1, dt=1ms
% load Noise10seriesT50sBL50Hz
% In=Signal;In(:,2:end)=Signal(:,2:end)*sc;

%-----------1 acao--------------------------------

% %% acao EC8 DNA Zona 1.1 Solo D Categ. II
% load EC8DNA11DIIacel.txt  %sinal sismo - valores em m/s2
% In(:,1)=EC8DNA11DIIacel(:,1);In(:,2)=EC8DNA11DIIacel(:,2)*sc;  

% %% acao EC8 DNA Zona 2.1 Solo D Categ. II
% load EC8DNA21DIIacel.txt  %sinal sismo - valores em m/s2
% In(:,1)=EC8DNA21DIIacel(:,1);In(:,2)=EC8DNA21DIIacel(:,2)*sc;

% %% acao EC8
% load SinalAcel.txt  %sinal sismo - valores em g
% In(:,1)=SinalAcel(:,1);In(:,2)=SinalAcel(:,2)*9.86*sc;

%% Sismo Elcentro
% cl=3; %2-FN;3-FP
load elcentro.txt
In=elcentro;
In(:,1)=elcentro(:,1);In(:,2)=elcentro(:,cl)*sc;

%%

% full_filename = fullfile(save_folder, filename+".txt");
% try
%     py_output = py.TXT_to_LTF.txt_to_ltf(full_filename_acq,save_folder);
%     % py_output is a Python string; convert to MATLAB char:
%     output_path = char(py_output);
%     fprintf('Python function returned output path: %s\n', output_path);
% catch ME
%     disp('Error calling Python function:');
%     disp(ME.message);
% end

%% Sismo Erzikan
% cl=3; %2-NS;3-EW
load erzikan.txt
In=erzikan;
In(:,1)=erzikan(:,1);In(:,2)=erzikan(:,cl)*sc;

%% Sismo Jiji
% cl=3; %2-NS;3-EW
load jiji.txt
In=jiji;
In(:,1)=jiji(:,1);In(:,2)=jiji(:,cl)*sc;

%% Sismo Kobe
% cl=3; %2-NS(FN);3-EW(FP)
load kobe.txt
In=kobe;
In(:,1)=kobe(:,1);In(:,2)=kobe(:,cl)*sc;

%% Sismo Newhall
% cl=3; %2-FN;3-FP
load newhall.txt
In=newhall;
In(:,1)=newhall(:,1);In(:,2)=newhall(:,cl)*sc;

%% Sismo Rinaldi
% cl=3; %2-FN;3-FP
load rinaldi.txt
In=rinaldi;
In(:,1)=rinaldi(:,1);In(:,2)=rinaldi(:,cl)*sc;

%% Sismo Sylmar
% cl=3; %2-FN;3-FP
load sylmar.txt
In=sylmar;
In(:,1)=sylmar(:,1);In(:,2)=sylmar(:,cl)*sc;
