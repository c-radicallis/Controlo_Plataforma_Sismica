%% Valores de Pico
clear all
close all

%todas as series
%
%file='TDOF_NoiseW1'; sensa=[10*0.6020*8 10/(4.1635/2) 10/(4.0192/2) -1 -1 1 -1]; %ganho 1
file='TDOF_NoiseW1_04';sensa=ones(7,1);   %ganho:09,09,07,06,05

%Escolher o inicio e no fim t(s)
ti=0;
tf='end'; 
%--------------------------------------------------------------------------
load(file);

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

xg=eval(['data.',sname{1},'.Y(4).Data']);xg=double(xg)*sensa(4);
xgr=eval(['data.',sname{1},'.Y(5).Data']);xgr=double(xgr)*sensa(5);
xig=eval(['data.',sname{1},'.Y(6).Data']);xig=double(xig)*sensa(6);
xsi=eval(['data.',sname{1},'.Y(7).Data']);xsi=double(xsi)*sensa(7);
ag=eval(['data.',sname{1},'.Y(1).Data']);ag=double(ag)*sensa(1);
ai=eval(['data.',sname{1},'.Y(2).Data']);ai=double(ai)*sensa(2);
as=eval(['data.',sname{1},'.Y(3).Data']);as=double(as)*sensa(3);



dimv=size(xg);
if dimv(1)==1
    time=time';xg=xg';xgr=xgr';xig=xig';xsi=xsi';ag=ag';ai=ai';as=as';
end
%
%iniciar em zero
    xg=xg-xg(1);xgr=xgr-xgr(1);xig=xig-xig(1);xsi=xsi-xsi(1);
%     ag=ag-mean(ag);ai=ai-mean(ai);as=as-mean(as);
    ag=ag-ag(1);ai=ai-ai(1);as=as-as(1);


%
if strcmp(tf,'end')
    clear tf
    tf=time(end);
end
    
nti=ceil(ti/dt)+1;
ntf=ceil(tf/dt)+1;


%----Modelo Experimental
timeexp=time(nti:ntf);
agexp=ag(nti:ntf);
aiexp=ai(nti:ntf);
asexp=as(nti:ntf);
xgexp=xg(nti:ntf);
xgrexp=xgr(nti:ntf);
xigexp=xig(nti:ntf);
xsiexp=xsi(nti:ntf);

xsgexp=xsiexp+xigexp;

%velocidades relativas
vgexp=diff(xgexp)/dt;
vigexp=diff(xigexp)/dt;
vsiexp=diff(xsiexp)/dt;
vsgexp=vsiexp+vigexp;

% figure('Position',[50 0 1000 800]),
% subplot(321),plot(timeexp'+ti,xgexp',timeexp+ti,xgrexp),ylabel('x_{gr} (mm)')
% subplot(323),plot(timeexp'+ti,xigexp'),ylabel('x_{ig} (mm)')
% subplot(325),plot(timeexp'+ti,xsiexp'),ylabel('x_{si} (mm)')
% subplot(322),plot(timeexp'+ti,agexp'),ylabel('a_{g} (m/s^2)')
% subplot(324),plot(timeexp'+ti,aiexp'),ylabel('a_{i} (m/s^2)'),xlabel('t (s)')
% subplot(326),plot(timeexp'+ti,asexp'),ylabel('a_{s} (m/s^2)'),xlabel('t (s)')


%Valores de pico
xgp=max(abs(xgexp));
vgp=max(abs(vgexp));
agp=max(abs(agexp));

out=[xgp/1000;vgp/1000;agp];open('out')

figure,
subplot(311),plot(timeexp'+ti,xgexp'),ylabel('x_{g} (mm)')
subplot(312),plot(timeexp(1:end-1)'+ti,vgexp'),ylabel('v_{g} (mm/s)')
subplot(313),plot(timeexp'+ti,agexp'),ylabel('a_{g} (m/s^2)')
