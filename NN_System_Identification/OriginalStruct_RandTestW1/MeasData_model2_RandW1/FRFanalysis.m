%% Visualizar dados

clear all
%close all

%-------ficheiros-----------------
%Ficheiros separados novos


file='TDOF_NoiseW1'; sensa=[10*0.6020*8 10/(4.1635/2) 10/(4.0192/2) -1 -1 1 -1]; %ganho 1
%file='TDOF_NoiseW1_09';sensa=ones(7,1);   %ganho:09,09,07,06,05


%sensa: Sensibilidade dos sensores:1-Ag;2-Ai;3-As;4-Xg;5-Xgref;6-Xig;7-Xsi;

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


%% FRF

close all

%Gravar FRF & Coherence
grab=1;


%close all

%----------------------
%TimeRange
tri=time(1);      %tempo de inicio
ttot=time(end);   %tempo total do sinal
%
%parametros para FRF
tr=[tri tri+ttot];
npt=10^4;%1.5*2^18;%ceil(ttot/dt); %nº de pontos
tw=1; %janela hanning
frange=[0.15 10]; %toda a gama [df 1/dt]
Hest=1; %estimador H1 para gravacao
%inde=0;  %integrar/derivar

%sensa=1./[4.2853 4.1635 4.0192]; %Inverso Sensibilidade Acelerometros
%sensa=1./[4.2853 3.7 4.0192];
sensa=1./[1 1 1];

[hXgrXg fo]=FreqResponseH3(time,xgr,xg,tr,npt,tw,frange,Hest,0,0,0);
%[hXgrAg fo]=FreqResponseH3(timefd',xgrfd',agfd',tr,npt,tw,frange,Hest,-2,0,0);

[hXgXig fo]=FreqResponseH3(time,xg,xig,tr,npt,tw,frange,Hest,-2,0,0);
[hXgXsi fo]=FreqResponseH3(time,xg,xsi,tr,npt,tw,frange,Hest,-2,0,0);
[hXgAg fo]=FreqResponseH3(time,xg,ag*sensa(1),tr,npt,tw,frange,Hest,-2,0,0);
[hAgXig fo]=FreqResponseH3(time,ag*sensa(1),xig,tr,npt,tw,frange,Hest,0,0,0);
[hAgXsi fo]=FreqResponseH3(time,ag*sensa(1),xsi,tr,npt,tw,frange,Hest,0,0,0);
[hAgAi fo]=FreqResponseH3(time,ag*sensa(1),ai*sensa(2),tr,npt,tw,frange,Hest,0,0,0);
[hAgAs fo]=FreqResponseH3(time,ag*sensa(1),as*sensa(3)',tr,npt,tw,frange,Hest,0,0,0);
[hXgXi fo]=FreqResponseH3(time,xg,xig+xg,tr,npt,tw,frange,Hest,0,0,0);
[hXgXs fo]=FreqResponseH3(time,xg,xsi+xig+xg,tr,npt,tw,frange,Hest,0,0,0);

[hXgrXgc fo]=FreqResponseH3(time,xgr,xg,tr,npt,tw,frange,0,0,0,0);
%[hXgrAgc fo]=FreqResponseH3(timefd',xgrfd',agfd',tr,npt,tw,frange,0,-2,0,0);

[hXgXigc fo]=FreqResponseH3(time,xg,xig,tr,npt,tw,frange,0,-2,0,0);
[hXgXsic fo]=FreqResponseH3(time,xg,xsi,tr,npt,tw,frange,0,-2,0,0);
[hXgAgc fo]=FreqResponseH3(time,xg,ag,tr,npt,tw,frange,0,-2,0,0);
[hAgXigc fo]=FreqResponseH3(time,ag,xig,tr,npt,tw,frange,0,0,0,0);
[hAgXsic fo]=FreqResponseH3(time,ag,xsi,tr,npt,tw,frange,0,0,0,0);
[hAgAic fo]=FreqResponseH3(time,ag,ai,tr,npt,tw,frange,0,0,0,0);
[hAgAsc fo]=FreqResponseH3(time,ag,as,tr,npt,tw,frange,0,0,0,0);
[hXgXic fo]=FreqResponseH3(time,xg,xig+xg,tr,npt,tw,frange,0,0,0,0);
[hXgXsc fo]=FreqResponseH3(time,xg,xsi+xig+xg,tr,npt,tw,frange,0,0,0,0);

%
%close all
figure('Position',[300 10 600 700]),
subplot(311),loglog(fo,abs(hXgrXg)),xlim(frange),ylabel('Mag'),title('Xg/X_{gr}')
subplot(312),semilogx(fo,angle(hXgrXg)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),semilogx(fo,abs(hXgrXgc)),xlim(frange),ylabel('Coherence'),xlabel('f (Hz)')


figure('Position',[300 10 600 700]),
subplot(311),loglog(fo,abs(hXgAg)),xlim(frange),ylabel('Mag'),title('Ag/Xg 1/s^2')
subplot(312),semilogx(fo,angle(hXgAg)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),semilogx(fo,abs(hXgAgc)),xlim(frange),ylabel('Coherence'),xlabel('f (Hz)')

figure('Position',[300 10 600 700]),
subplot(311),loglog(fo,abs(hXgXig),fo,abs(hAgXig)/1000),xlim(frange),title('Xig/Ag'),ylabel('Mag (m/m/s^2)')
subplot(312),semilogx(fo,angle(hXgXig)*180/pi,fo,angle(hAgXig)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),semilogx(fo,abs(hXgXigc),fo,abs(hAgXigc)),xlim(frange),ylabel('Coherence')
xlabel('f (Hz)'),legend('Xig/Xg 1/s^2','Xig/Ag')


figure('Position',[300 10 600 700]),
subplot(311),loglog(fo,abs(hXgXsi),fo,abs(hAgXsi)/1000),xlim(frange),title('Xsi/Ag'),ylabel('Mag (m/m/s^2)')
subplot(312),semilogx(fo,angle(hXgXsi)*180/pi,fo,angle(hAgXsi)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),semilogx(fo,abs(hXgXsic),fo,abs(hAgXsic)),xlim(frange),ylabel('Coherence')
xlabel('f (Hz)'),legend('Xsi/Xg 1/s^2','Xsi/Ag')

figure('Position',[300 10 600 700]),
subplot(311),loglog(fo,abs(hXgXi),fo,abs(hAgAi)),xlim(frange),title('Xi/Xg'),ylabel('Mag')
subplot(312),semilogx(fo,angle(hXgXi)*180/pi,fo,angle(hAgAi)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),semilogx(fo,abs(hXgXic),fo,abs(hAgAic)),xlim(frange),ylabel('Coherence')
xlabel('f (Hz)'),legend('(Xig+Xg)/Xg','Ai/Ag')

figure('Position',[300 10 600 700]),
subplot(311),loglog(fo,abs(hXgXs),fo,abs(hAgAs)),xlim(frange),title('Xs/Xg'),ylabel('Mag')
subplot(312),semilogx(fo,angle(hXgXs)*180/pi,fo,angle(hAgAs)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),semilogx(fo,abs(hXgXsc),fo,abs(hAgAsc)),xlim(frange),ylabel('Coherence')
xlabel('f (Hz)'),legend('(Xsi+Xig+Xg)/Xg','As/Ag')



% Gravar FRF

if grab==1
    
FRF.f=fo;
FRF.AgS2Xg=hXgAg;
FRF.XigS2Xg=hXgXig;
FRF.XsiS2Xg=hXgXsi;
FRF.XiXg=hXgXi;
FRF.XsXg=hXgXs;
FRF.XigAg=hAgXig;
FRF.XsiAg=hAgXsi;
FRF.AiAg=hAgAi;
FRF.AsAg=hAgAs;

CO.AgS2Xg=hXgAgc;
CO.XigS2Xg=hXgXigc;
CO.XsiS2Xg=hXgXsic;
CO.XiXg=hXgXic;
CO.XsXg=hXgXsc;
CO.XigAg=hAgXigc;
CO.XsiAg=hAgXsic;
CO.AiAg=hAgAic;
CO.AsAg=hAgAsc;

uisave({'FRF','CO'},'2DOF_FreqResp_G')
end

%% Espetro
%
%Sinal
timesig=time;
ysig=as;

%FFT
nptos=length(timesig);%2*10^4;   %nº pontos FFT
ov=50;
tw=1;


%FFT
%uma serie
SignalSpectra([timesig ysig],nptos,ov,tw,2)


%varias series
nptost=nptos;%nptost=ceil(nptos*[1,1,1/resr]);
%sinal(1).data=[timesig,ysig];
%sinal(2).data=[timesf,ysigf];

SignalSpectraComp(sinal,nptost,ov,tw,2)


%integrar o sinal
%Integra([timesig ysig],0,0)