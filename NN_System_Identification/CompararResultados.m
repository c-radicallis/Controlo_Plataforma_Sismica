%Identificação via dominio do tempo
%% Visualizar dados

clear all
%close all

%-------ficheiros-----------------
%Ficheiros separados novos

%file='TDOF_SigOut';
%file='TDOF_SigOutFD';

file='TDOF_PassFD_Sylmar_s2_100_V5.mat';

%
%-------dados-----------------

data=load(file);
sname=fieldnames(data);

%
time=eval(['data.',sname{1},'.X.Data']);
dt=mean(diff(time));
%
%saida '*.Y(i).Data':
%i:1-Ag;2-Ai;3-As;4-Xg;5-Xgref;6-Xig;7-Xsi;

xg=eval(['data.',sname{1},'.Y(4).Data']);xg=double(xg);
xgr=eval(['data.',sname{1},'.Y(5).Data']);xgr=double(xgr);
xig=eval(['data.',sname{1},'.Y(6).Data']);xig=double(xig);
xsi=eval(['data.',sname{1},'.Y(7).Data']);xsi=double(xsi);
ag=eval(['data.',sname{1},'.Y(1).Data']);ag=double(ag);
ai=eval(['data.',sname{1},'.Y(2).Data']);ai=double(ai);
as=eval(['data.',sname{1},'.Y(3).Data']);as=double(as);

dimv=size(xg);
if dimv(1)==1
    time=time';xg=xg';xgr=xgr';xig=xig';xsi=xsi';ag=ag';ai=ai';as=as';
end
%
figure('Position',[50 0 1000 800]),
subplot(321),plot(time,xg,time,xgr),ylabel('x_{gr} (mm)')
subplot(323),plot(time,xig),ylabel('x_{ig} (mm)')
subplot(325),plot(time,xsi),ylabel('x_{si} (mm)')
subplot(322),plot(time,ag),ylabel('a_{g} (m/s^2)')
subplot(324),plot(time,ai),ylabel('a_{i} (m/s^2)'),xlabel('t (s)')
subplot(326),plot(time,as),ylabel('a_{s} (m/s^2)'),xlabel('t (s)')


%% filtrar e decimar

clear dataf

%input:Data;Filteer;Fdecimation;Inic0(0,1);Mostrar(0,1,2);Gravar(0,1);
FilterData(data,'FilterLP10',100,1,2,1)



%% Comparacao Encostado/desencostado
%Comparacao de Sinais

clear all
%close all

%iniciar em zero
inic=2;  %0-nao altera inicio;1-iniciar em zero;2-media em 0;

%diferenca de tempo


%-------ficheiros-----------------
%saida '*.Y(i).Data':
%i:1-Ag;2-Ai;3-As;4-Xg;5-Xgref;6-Xig;7-Xsi;

%Rand Tests c/massa adicional 300kg na base
%input: RandF02to100hzT100s

%ficheiros
% file{1}='TDOF_SigOut';   %sem filtragem/decimacao
% file{2}='TDOF_SigOutFilDec';  %Filt&Dec200hz pos-processament
% file{3}='TDOF_SigOutFilDec100hz';  %Filt&Dec100hz pos-processament
% file{4}='TDOF_SigOutFD';  %Filt&Dec online

file{1}='TDOF_SigOut_A1s1';   %sem filtragem/decimacao
file{2}='TDOF_SigOut_A1s1FilDec';  %Filt&Dec200hz pos-processament
file{3}='TDOF_SigOutFD_A1s1';  %Filt&Dec online


dtime=[0 -0.132 -0.132];

figure('Position',[50 0 1000 800]),
for i=1:3

clear data sname time xg xgr xig xsi ag ai as

%-------dados-----------------
data=load(file{i});
sname=fieldnames(data);

%
time=eval(['data.',sname{1},'.X.Data']);
dt=mean(diff(time));
time=time+dtime(i);


%saida '*.Y(i).Data':
%i:1-Ag;2-Ai;3-As;4-Xg;5-Xgref;6-Xig;7-Xsi;
xg=eval(['data.',sname{1},'.Y(4).Data']);xg=double(xg);
%xgr=eval(['data.',sname{1},'.Y(5).Data']);xgr=double(xgr);
xig=eval(['data.',sname{1},'.Y(6).Data']);xig=double(xig);
xsi=eval(['data.',sname{1},'.Y(7).Data']);xsi=double(xsi);
ag=eval(['data.',sname{1},'.Y(1).Data']);ag=double(ag);
ai=eval(['data.',sname{1},'.Y(2).Data']);ai=double(ai);
as=eval(['data.',sname{1},'.Y(3).Data']);as=double(as);

if inic==1
    xg=xg-xg(1);%xgr=xgr-xgr(1);
    xig=xig-xig(1);
    xsi=xsi-xsi(1);
    ag=ag-ag(1);
    ai=ai-ai(1);
    as=as-as(1);
elseif inic==2
    xg=xg-mean(xg);%xgr=xgr-xgr(1);
    xig=xig-mean(xig);
    xsi=xsi-mean(xsi);
    ag=ag-mean(ag);
    ai=ai-mean(ai);
    as=as-mean(as);
end

subplot(321),
hold all
plot(time',xg'),ylabel('x_{gr} (mm)')
hold off

subplot(322),
hold all
plot(time',ag'),ylabel('a_{g} (m/s^2)')
hold off

subplot(323),
hold all
plot(time',xig'),ylabel('x_{ig} (mm)')
hold off

subplot(324),
hold all
plot(time',ai'),ylabel('a_{i} (m/s^2)'),xlabel('t (s)')
hold off

subplot(325),
hold all
plot(time',xsi'),ylabel('x_{si} (mm)')
hold off

subplot(326),
hold all
plot(time',as'),ylabel('a_{s} (m/s^2)'),xlabel('t (s)')
hold off

end
