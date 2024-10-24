%Ruido Filtrado
clear all, close all
s=tf('s');
%funcao para gerar series de ruido filtrado
%especificando as entradas, projeta-se um filtro de ordem 'n' que depois
%sera usado para filtrar uma serie de ruido branco no simulink; na serie
%final sao aplicadas as funsoes cosine-taper nas extremidas; o processo
%e iterativo; depois de observada a serie e o seu espetro, e sendo esses os
%pretendidos, grava-se a serie;

%entradas:-----------------------------------------------------------------
fn1f=0.1;   %fn1f:frequencia do passa-alto
fn2f=30;    %fn2f:frequencia do passa-baixo
pf=3;       %pf:peso nas freq. fn1f e fn2f - f1=fn1f/pf;f2=fn2f*pf; 
fc=2;       %fc:frequencia central-freq. intermedia para ajuste do ganho do filtro
od=16;      %od:ordem do filtro
dur=20;     %dur:duracao do sinal
dt=10^-3;   %dt:discretizacao no tempo
gain=0.01;  %gain:ganho a aplicar ao sinal no tempo
grava=0;    %gravar a serie:grava=0-nao; grava=1-sim
%--------------------------------------------------------------------------

%carateristicas da serie
tsub=2;     %tempo de transição para as funcoes cosine-taper
Intime=0;   %tempo de inicio da simulacao
Sttime=dur; %tempo final da simulacao

zetaf=0.77;  %amortecimento do filtro
wn1f=2*pi*fn1f/pf; %freq em rad/s
wn2f=2*pi*fn2f*pf;  

%funcao de transferencia-filtro
FTfp=minreal(s^od/(s^2+2*zetaf*wn1f*s+wn1f^2)^(od/2)/(s^2+2*zetaf*wn2f*s+wn2f^2)^(od/2));
%FTfp=minreal(s^2/(s^2+2*zetaf*wn1f*s+wn1f^2)/(s^2+2*zetaf*wn2f*s+wn2f^2));

%ganho na frequencia intermedia entre fn1 e fn2
[mag,phase] = bode(FTfp,fc*2*pi);
FTf=1/mag*FTfp; %ajuste do ganho

%bode diagram
P = bodeoptions; 
P.PhaseVisible = 'off';
P.FreqUnits = 'Hz';
P.Grid='on';
P.Title.String='Filter';
P.XLim={[fn1f/10 fn2f*10]};
%figure,bode(FTf,P)

%simulacao: serie
options = simset('FixedStep',dt,'Solver','ode3'); %p. integrcao fixo
sim('RuidoFiltrado',[Intime Sttime],options)

tv=rf.time;
val=rf.signals.values;

%valor RMS
rmsval=sqrt(mean(val.^2));

%Estimador com janela Tukey (cosine taper) com r=0.05
le=8192;        %nº de pontos de cada segmento de sinal para estimar a resposta em frequencia
ov=30;          %percentagem de sobreposicao dos segmentos
%para obter o periodograma só com a janela escolher:
%le=nº total de pontos do sinal; ov=0;  

[msf,msd]=espectrowelch([tv val],le,ov,0);

%determinacao da amp media no intervalo de ferquencias escolhido
df=msf(end)-msf(end-1); %discretizacao em frequencia
n1=ceil(fn1f/df);
n2=ceil(fn2f/df);
medaf=mean(msd(n1:n2)); %media

if grava==0 %nao mostra o grfico qd gravar
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)]),
subplot(2,2,[1 3]),bode(FTf,P)
subplot(2,2,2),plot(tv,val,tv,rmsval*ones(size(val))),
title('time signal'),xlabel('t (s)'),ylabel('A')
legend('signal', ['RMS value =',sprintf('%5.4f',rmsval)])
subplot(2,2,4),semilogx(msf,msd,msf,medaf*ones(size(msd))),XLim([fn1f/pf/2 fn2f*pf*2])
title('Welch Spectrum'),xlabel('f (Hz)'),ylabel('|A|'),
legend('spectrum', ['mean value =',sprintf('%5.4f',medaf)])

else
% Gravar o ficheiro
Ruido=[tv val];
uisave('Ruido','RuidoFiltrado')
end
%% Comparacao de filtros

clear all, close all
s=tf('s');
fn1f=0.1;   %fn1f:frequencia do passa-alto
fn2f=30;    %fn2f:frequencia do passa-baixo
pf=3;       %pf:peso nas freq. fn1f e fn2f - f1=fn1f/pf;f2=fn2f*pf; 
fc=2;       %fc:frequencia central-freq. intermedia para ajuste do ganho do filtro
od=16;      %od:ordem do filtro
zetaf=0.77;  %amortecimento do filtro
wn1f=2*pi*fn1f/pf; %freq em rad/s
wn2f=2*pi*fn2f*pf;  

FTfpn=minreal(s^od/(s^2+2*zetaf*wn1f*s+wn1f^2)^(od/2)/(s^2+2*zetaf*wn2f*s+wn2f^2)^(od/2));
[mag,phase] = bode(FTfpn,fc*2*pi);
FTfn=1/mag*FTfpn; %ajuste do ganho

FTfp2=minreal(s^2/(s^2+2*zetaf*wn1f*s+wn1f^2)/(s^2+2*zetaf*wn2f*s+wn2f^2));
[mag2,phase2] = bode(FTfp2,fc*2*pi);
FTf2=1/mag2*FTfp2; %ajuste do ganho

%bode diagram
P = bodeoptions; 
P.PhaseVisible = 'off';
P.FreqUnits = 'Hz';
P.Grid='on';
P.Title.String='Filter';
P.XLim={[fn1f/10 fn2f*10]};
figure,bode(FTfn,FTf2,P)
