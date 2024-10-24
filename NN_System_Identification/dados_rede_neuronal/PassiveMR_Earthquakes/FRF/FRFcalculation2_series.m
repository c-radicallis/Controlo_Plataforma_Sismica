%% Visualizar dados

clear all
close all

%-------ficheiros-----------------

gain='Gdif'; %Ganho: '100'-1;'90'-0.9;...;'50'-0.5;'Gdiff'-ganhos diferentes
quake='A1&2'; %Acao: 'A1';'A2';%A1&2
Vin='5';  %tensao eletrica:0-5
Serie='all'; %serie:1-10; all;
%
%------------------------------------------------------------
file=strcat('TDOF_PassFD_',quake,'_',Serie,'_',gain,'_V',Vin);

sensa=ones(8,1);  %acerto dos sinais


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

xg=eval(['data.',sname{1},'.Y(5).Data']);xg=double(xg)*sensa(5);
xgr=eval(['data.',sname{1},'.Y(6).Data']);xgr=double(xgr)*sensa(6);
xig=eval(['data.',sname{1},'.Y(7).Data']);xig=double(xig)*sensa(7);
xsi=eval(['data.',sname{1},'.Y(8).Data']);xsi=double(xsi)*sensa(8);
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
npt=4096;%1.5*2^18;%ceil(ttot/dt); %nº de pontos
%npt=10^4;
frange=[0.15 5]; %toda a gama [df 1/dt]
Hest=1; %estimador H1 para gravacao

%janela de dados
%wd=tukeywin(npt,1);  %{'tukey',1}={hanning},{'tukey',0}={rectangular}
%wd=rectwin(npt);
wd=hanning(npt);
%wd=hamming(npt);

%overlap
ov=2/3;  %1/2 ou 2/3  

esc='lin'; %escala

%sensa=1./[4.2853 4.1635 4.0192]; %Inverso Sensibilidade Acelerometros
%sensa=1./[4.2853 3.7 4.0192];
sensa=1./[1 1 1];

[hXgrXg fo]=FreqResponseH4(time,xgr,xg,tr,wd,ov,frange,Hest,0,0,0);
%[hXgrAg fo]=FreqResponseH4(timefd',xgrfd',agfd',tr,npt,tw,frange,Hest,-2,0,0);

[hXgXig fo]=FreqResponseH4(time,xg,xig,tr,wd,ov,frange,Hest,-2,0,0);
[hXgXsi fo]=FreqResponseH4(time,xg,xsi,tr,wd,ov,frange,Hest,-2,0,0);
[hXgAg fo]=FreqResponseH4(time,xg,ag*sensa(1),tr,wd,ov,frange,Hest,-2,0,0);
[hAgXig fo]=FreqResponseH4(time,ag*sensa(1),xig,tr,wd,ov,frange,Hest,0,0,0);
[hAgXsi fo]=FreqResponseH4(time,ag*sensa(1),xsi,tr,wd,ov,frange,Hest,0,0,0);
[hAgAi fo]=FreqResponseH4(time,ag*sensa(1),ai*sensa(2),tr,wd,ov,frange,Hest,0,0,0);
[hAgAs fo]=FreqResponseH4(time,ag*sensa(1),as*sensa(3)',tr,wd,ov,frange,Hest,0,0,0);
[hXgXi fo]=FreqResponseH4(time,xg,xig+xg,tr,wd,ov,frange,Hest,0,0,0);
[hXgXs fo]=FreqResponseH4(time,xg,xsi+xig+xg,tr,wd,ov,frange,Hest,0,0,0);

[hXgrXgc fo]=FreqResponseH4(time,xgr,xg,tr,wd,ov,frange,0,0,0,0);
%[hXgrAgc fo]=FreqResponseH4(timefd',xgrfd',agfd',tr,wd,ov,frange,0,-2,0,0);

[hXgXigc fo]=FreqResponseH4(time,xg,xig,tr,wd,ov,frange,0,-2,0,0);
[hXgXsic fo]=FreqResponseH4(time,xg,xsi,tr,wd,ov,frange,0,-2,0,0);
[hXgAgc fo]=FreqResponseH4(time,xg,ag,tr,wd,ov,frange,0,-2,0,0);
[hAgXigc fo]=FreqResponseH4(time,ag,xig,tr,wd,ov,frange,0,0,0,0);
[hAgXsic fo]=FreqResponseH4(time,ag,xsi,tr,wd,ov,frange,0,0,0,0);
[hAgAic fo]=FreqResponseH4(time,ag,ai,tr,wd,ov,frange,0,0,0,0);
[hAgAsc fo]=FreqResponseH4(time,ag,as,tr,wd,ov,frange,0,0,0,0);
[hXgXic fo]=FreqResponseH4(time,xg,xig+xg,tr,wd,ov,frange,0,0,0,0);
[hXgXsc fo]=FreqResponseH4(time,xg,xsi+xig+xg,tr,wd,ov,frange,0,0,0,0);

%
%close all
if strcmp(esc,'log')

figure('Position',[300 10 600 700]),
subplot(311),
loglog(fo,abs(hXgrXg)),xlim(frange),ylabel('Mag'),title('Xg/X_{gr}')
subplot(312),semilogx(fo,angle(hXgrXg)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),semilogx(fo,abs(hXgrXgc)),xlim(frange),ylabel('Coherence'),xlabel('f (Hz)'),ylim([0 1])


figure('Position',[300 10 600 700]),
subplot(311),loglog(fo,abs(hXgAg)),xlim(frange),ylabel('Mag'),title('Ag/Xg 1/s^2')
subplot(312),semilogx(fo,angle(hXgAg)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),semilogx(fo,abs(hXgAgc)),xlim(frange),ylabel('Coherence'),xlabel('f (Hz)'),ylim([0 1])

figure('Position',[300 10 600 700]),
subplot(311),loglog(fo,abs(hXgXig),fo,abs(hAgXig)/1000),xlim(frange),title('Xig/Ag'),ylabel('Mag (m/m/s^2)')
subplot(312),semilogx(fo,angle(hXgXig)*180/pi,fo,angle(hAgXig)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),semilogx(fo,abs(hXgXigc),fo,abs(hAgXigc)),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)'),legend('Xig/Xg 1/s^2','Xig/Ag')


figure('Position',[300 10 600 700]),
subplot(311),loglog(fo,abs(hXgXsi),fo,abs(hAgXsi)/1000),xlim(frange),title('Xsi/Ag'),ylabel('Mag (m/m/s^2)')
subplot(312),semilogx(fo,angle(hXgXsi)*180/pi,fo,angle(hAgXsi)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),semilogx(fo,abs(hXgXsic),fo,abs(hAgXsic)),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)'),legend('Xsi/Xg 1/s^2','Xsi/Ag')

figure('Position',[300 10 600 700]),
subplot(311),loglog(fo,abs(hXgXi),fo,abs(hAgAi)),xlim(frange),title('Xi/Xg'),ylabel('Mag')
subplot(312),semilogx(fo,angle(hXgXi)*180/pi,fo,angle(hAgAi)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),semilogx(fo,abs(hXgXic),fo,abs(hAgAic)),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)'),legend('(Xig+Xg)/Xg','Ai/Ag')

figure('Position',[300 10 600 700]),
subplot(311),loglog(fo,abs(hXgXs),fo,abs(hAgAs)),xlim(frange),title('Xs/Xg'),ylabel('Mag')
subplot(312),semilogx(fo,angle(hXgXs)*180/pi,fo,angle(hAgAs)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),semilogx(fo,abs(hXgXsc),fo,abs(hAgAsc)),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)'),legend('(Xsi+Xig+Xg)/Xg','As/Ag')

elseif strcmp(esc,'lin')
    
figure('Position',[300 10 600 700]),
subplot(311),
semilogy(fo,abs(hXgrXg)),xlim(frange),ylabel('Mag'),title('Xg/X_{gr}')
subplot(312),plot(fo,angle(hXgrXg)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),plot(fo,abs(hXgrXgc)),xlim(frange),ylabel('Coherence'),xlabel('f (Hz)'),ylim([0 1])


figure('Position',[300 10 600 700]),
subplot(311),semilogy(fo,abs(hXgAg)),xlim(frange),ylabel('Mag'),title('Ag/Xg 1/s^2')
subplot(312),plot(fo,angle(hXgAg)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),plot(fo,abs(hXgAgc)),xlim(frange),ylabel('Coherence'),xlabel('f (Hz)'),ylim([0 1])

figure('Position',[300 10 600 700]),
subplot(311),semilogy(fo,abs(hXgXig),fo,abs(hAgXig)/1000),xlim(frange),title('Xig/Ag'),ylabel('Mag (m/m/s^2)')
subplot(312),plot(fo,angle(hXgXig)*180/pi,fo,angle(hAgXig)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),plot(fo,abs(hXgXigc),fo,abs(hAgXigc)),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)'),legend('Xig/Xg 1/s^2','Xig/Ag')


figure('Position',[300 10 600 700]),
subplot(311),semilogy(fo,abs(hXgXsi),fo,abs(hAgXsi)/1000),xlim(frange),title('Xsi/Ag'),ylabel('Mag (m/m/s^2)')
subplot(312),plot(fo,angle(hXgXsi)*180/pi,fo,angle(hAgXsi)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),plot(fo,abs(hXgXsic),fo,abs(hAgXsic)),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)'),legend('Xsi/Xg 1/s^2','Xsi/Ag')

figure('Position',[300 10 600 700]),
subplot(311),semilogy(fo,abs(hXgXi),fo,abs(hAgAi)),xlim(frange),title('Xi/Xg'),ylabel('Mag')
subplot(312),plot(fo,angle(hXgXi)*180/pi,fo,angle(hAgAi)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),plot(fo,abs(hXgXic),fo,abs(hAgAic)),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)'),legend('(Xig+Xg)/Xg','Ai/Ag')

figure('Position',[300 10 600 700]),
subplot(311),semilogy(fo,abs(hXgXs),fo,abs(hAgAs)),xlim(frange),title('Xs/Xg'),ylabel('Mag')
subplot(312),plot(fo,angle(hXgXs)*180/pi,fo,angle(hAgAs)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),plot(fo,abs(hXgXsc),fo,abs(hAgAsc)),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)'),legend('(Xsi+Xig+Xg)/Xg','As/Ag')

end


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




%% Juntar varias series
%%

clear all
%close all

%Filtro passa-baixo
filter='none';  %'FilterLP10';'FilterLP30';'FilterLP40';


%Resample (2xfn)
fnr=0;

%Remover no inicio e no fim t(s)
trem=0;

%Inverso Sensibilidade Acelerometros
sensa=[1 1 1];

%-------ficheiros-----------------
%saida '*.Y(i).Data':
%i:1-Ag;2-Ai;3-As;4-Xg;5-Xgref;6-Xig;7-Xsi;

%serie=1:10;
%quake='A1'; gain={'100','100','100','90','90','80','90','80','90','80'};%'Acao1'
%quake='A2'; gain={'80','100','100','90','90','70','80','90','100','60'};%'Acao2'
%Vin='5';

serie=1:2;
Vin='5';
%
xgo=zeros(110*200*length(serie),1);
xgro=xgo;xigo=xgo;xsio=xgo;ago=xgo;aio=xgo;aso=xgo;
ind=1; %indice para colocar nova serie
for i=1:length(serie)

clear data sname timei xgi xgri xigi xsii agi aii asi time xg xgr xig xsi ag ai as ...
      timef XP xgf xgrf xigf xsif agf aif asf timefd xgfd xgrfd xigfd xsifd agfd aifd asfd  

%file=strcat('TDOF_PassFD_',quake,'_s',num2str(serie(i)),'_',gain{i},'_V',Vin);
  
file=strcat('TDOF_PassFD_A',num2str(serie(i)),'_all_Gdif_V',Vin);
%
%-------dados-----------------

data=load(file);
sname=fieldnames(data);
%
timei=eval(['data.',sname{1},'.X.Data']);
dt=mean(diff(timei));
%
%saida '*.Y(i).Data':
%i:1-Ag;2-Ai;3-As;5-Xg;6-Xgref;7-Xig;8-Xsi;
fai=eval(['data.',sname{1},'.Y(4).Data']);fai=double(fai);
xgi=eval(['data.',sname{1},'.Y(5).Data']);xgi=double(xgi);
xgri=eval(['data.',sname{1},'.Y(6).Data']);xgri=double(xgri);
xigi=eval(['data.',sname{1},'.Y(7).Data']);xigi=double(xigi);
xsii=eval(['data.',sname{1},'.Y(8).Data']);xsii=double(xsii);
agi=eval(['data.',sname{1},'.Y(1).Data']);agi=double(agi)*sensa(1);
aii=eval(['data.',sname{1},'.Y(2).Data']);aii=double(aii)*sensa(2);
asi=eval(['data.',sname{1},'.Y(3).Data']);asi=double(asi)*sensa(3);

%remover inicio e fim
nrem=ceil(trem/dt)+1;
time=timei(1:end-2*nrem+1);
fa=fai(nrem:end-nrem);
xg=xgi(nrem:end-nrem);
xgr=xgri(nrem:end-nrem);
xig=xigi(nrem:end-nrem);
xsi=xsii(nrem:end-nrem);
ag=agi(nrem:end-nrem);
ai=aii(nrem:end-nrem);
as=asi(nrem:end-nrem);

%iniciar em zero
% xg=xg-xg(1);
% xgr=xgr-xgr(1);
% xig=xig-xig(1);
% xsi=xsi-xsi(1);
% ag=ag-ag(1);
% ai=ai-ai(1);
% as=as-as(1);

%iniciar no valor medio
% xg=xg-mean(xg);
% xgr=xgr-mean(xgr);
% xig=xig-mean(xig);
% xsi=xsi-mean(xsi);
% ag=ag-mean(ag);
% ai=ai-mean(ai);
% as=as-mean(as);

%iniciar nem zero
xg=xg-xg(1);
xgr=xgr-xgr(1);
xig=xig-xig(1);
xsi=xsi-xsi(1);
ag=ag-ag(1);
ai=ai-ai(1);
as=as-as(1);

if strcmp(filter,'none')
    faf=fa';
    xgf=xg';
    xgrf=xgr';
    xigf=xig';
    xsif=xsi';
    agf=ag';
    aif=ai';
    asf=as';
    timef=time';
else
options = simset('FixedStep',dt,'Solver','FixedStepDiscrete'); %p. integrcao fixo
[timef, XP, xgf] = sim(filter, [time(1) time(end)], options,[time' xg']);
[timef, XP, xgrf] = sim(filter, [time(1) time(end)], options,[time' xgr']);
[timef, XP, xigf] = sim(filter, [time(1) time(end)], options,[time' xig']);
[timef, XP, xsif] = sim(filter, [time(1) time(end)], options,[time' xsi']);
[timef, XP, agf] = sim(filter, [time(1) time(end)], options,[time' ag']);
[timef, XP, aif] = sim(filter, [time(1) time(end)], options,[time' ai']);
[timef, XP, asf] = sim(filter, [time(1) time(end)], options,[time' as']);
end

if fnr==0
    fafd=faf;
    xgfd=xgf;
    xgrfd=xgrf;
    xigfd=xigf;
    xsifd=xsif;
    agfd=agf;
    aifd=aif;
    asfd=asf;
    timefd=timef;
    fnr2=1/dt;
else
resr=ceil(1/dt/fnr); %resample rate
xgfd = resample(xgf,1,resr);
xgrfd = resample(xgrf,1,resr);
xigfd = resample(xigf,1,resr);
xsifd = resample(xsif,1,resr);
agfd = resample(agf,1,resr);
aifd = resample(aif,1,resr);
asfd = resample(asf,1,resr);
timefd=(0:1:length(agfd)-1)'*1/fnr;
fnr2=fnr;
end

fao(ind:ind+length(timefd)-1)=fafd;
xgo(ind:ind+length(timefd)-1)=xgfd;
xgro(ind:ind+length(timefd)-1)=xgrfd;
xigo(ind:ind+length(timefd)-1)=xigfd;
xsio(ind:ind+length(timefd)-1)=xsifd;
ago(ind:ind+length(timefd)-1)=agfd;
aio(ind:ind+length(timefd)-1)=aifd;
aso(ind:ind+length(timefd)-1)=asfd;

ind=ind+length(timefd);
end


%output
out.X.Data=(0:1:ind-2)'/fnr2;
out.Y(1).Data=ago(1:ind-1);
out.Y(2).Data=aio(1:ind-1);
out.Y(3).Data=aso(1:ind-1);
out.Y(4).Data=fao(1:ind-1);
out.Y(5).Data=xgo(1:ind-1);
out.Y(6).Data=xgro(1:ind-1);
out.Y(7).Data=xigo(1:ind-1);
out.Y(8).Data=xsio(1:ind-1);


%
figure('Position',[50 0 1000 800]),
subplot(321),plot(out.X.Data,out.Y(5).Data,out.X.Data,out.Y(6).Data),ylabel('x_{gr} (mm)')
subplot(323),plot(out.X.Data,out.Y(7).Data),ylabel('x_{ig} (mm)')
subplot(325),plot(out.X.Data,out.Y(8).Data),ylabel('x_{si} (mm)')
subplot(322),plot(out.X.Data,out.Y(1).Data),ylabel('a_{g} (m/s^2)')
subplot(324),plot(out.X.Data,out.Y(2).Data),ylabel('a_{i} (m/s^2)'),xlabel('t (s)')
subplot(326),plot(out.X.Data,out.Y(3).Data),ylabel('a_{s} (m/s^2)'),xlabel('t (s)')


uisave('out',strcat('TDOF_PassFD_A_all_Gdif_V'))

%% Valores de Pico
%%
clear all
close all

%-------ficheiros-----------------

gain='1'; %'1'-1;'_09'-0.9;...;'_05'-0.5;
quake='A2'; %'A1';'A2';
%
%------------------------------------------------------------
file=strcat('TDOF_SigOutFD_',quake,'_all_',gain);
sensa=[1 1 1 -1 -1 1 -1];
%sensa=ones(7,1);

%sensa: Sensibilidade dos sensores:1-Ag;2-Ai;3-As;4-Xg;5-Xgref;6-Xig;7-Xsi;

%Escolher o inicio e no fim t(s)
ti=0;
tf='end'; 
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
ntf=ceil(tf/dt);


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
