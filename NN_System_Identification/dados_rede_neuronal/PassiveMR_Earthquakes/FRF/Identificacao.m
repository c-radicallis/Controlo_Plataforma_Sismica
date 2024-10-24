% Identificacao
%%

clear all 
close all

%ficheiro com FRFs
gain='diff'; %'1'-1;'_09'-0.9;...;'_05'-0.5; 'diff'-ganhos diferentes
quake='A1&2'; %'A1';'A2';
npt='4096';     %n� de pontos
Vin='2';

%--------------------------------------------------------------------
data=load(strcat('TDOF_FreqResp_',quake,'G',gain,'V',Vin,'_',npt));

%gama de frequencias
frange=[0.2 5];
esc='lin';

%Modelo espacial
%Massa/Rigidez
mk='m';  %'m'-massa conhecida; 'k'-rigidez conhecida

%Amortecimento: 1) alpha & beta ~=0; 2) beta=0;
y0 = [1 1e-1];lby=y0*10^-3;uby=y0*10^3;                 
%y0 = [0.1 0];lby=y0*10^-3;uby=y0*10^3; Aeqop=[0 1];beqop=0; %c/beta=0;

%--------------------------------------------------------------------------
   
figure('Position',[300 10 600 700]),
subplot(311),
semilogy(data.FRF.f,abs(data.FRF.AgS2Xg)),xlim(frange),ylabel('Mag'),title('Xg/X_{gr}')
subplot(312),plot(data.FRF.f,angle(data.FRF.AgS2Xg)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (�)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),plot(data.FRF.f,data.CO.AgS2Xg),xlim(frange),ylabel('Coherence'),xlabel('f (Hz)'),ylim([0 1])


figure('Position',[300 10 600 700]),
subplot(311),semilogy(data.FRF.f,abs(data.FRF.XigS2Xg),data.FRF.f,abs(data.FRF.XigAg/1000)),xlim(frange),title('Xig/Ag'),ylabel('Mag (m/m/s^2)')
subplot(312),plot(data.FRF.f,angle(data.FRF.XigS2Xg)*180/pi,data.FRF.f,angle(data.FRF.XigAg)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (�)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),plot(data.FRF.f,data.CO.XigS2Xg,data.FRF.f,data.CO.XigAg),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)'),legend('Xig/Xg 1/s^2','Xig/Ag')


figure('Position',[300 10 600 700]),
subplot(311),semilogy(data.FRF.f,abs(data.FRF.XsiS2Xg),data.FRF.f,abs(data.FRF.XsiAg/1000)),xlim(frange),title('Xsi/Ag'),ylabel('Mag (m/m/s^2)')
subplot(312),plot(data.FRF.f,angle(data.FRF.XsiS2Xg)*180/pi,data.FRF.f,angle(data.FRF.XsiAg)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (�)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),plot(data.FRF.f,data.CO.XsiS2Xg,data.FRF.f,data.CO.XsiAg),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)'),legend('Xsi/Xg 1/s^2','Xsi/Ag')

figure('Position',[300 10 600 700]),
subplot(311),semilogy(data.FRF.f,abs(data.FRF.XiXg),data.FRF.f,abs(data.FRF.AiAg)),xlim(frange),title('Xi/Xg'),ylabel('Mag')
subplot(312),plot(data.FRF.f,angle(data.FRF.XiXg)*180/pi,data.FRF.f,angle(data.FRF.AiAg)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (�)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),plot(data.FRF.f,data.CO.XiXg,data.FRF.f,data.CO.AiAg),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)'),legend('(Xig+Xg)/Xg','Ai/Ag')

figure('Position',[300 10 600 700]),
subplot(311),semilogy(data.FRF.f,abs(data.FRF.XsXg),data.FRF.f,abs(data.FRF.AsAg)),xlim(frange),title('Xs/Xg'),ylabel('Mag')
subplot(312),plot(data.FRF.f,angle(data.FRF.XsXg)*180/pi,data.FRF.f,angle(data.FRF.AsAg)*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (�)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
subplot(313),plot(data.FRF.f,data.CO.XsXg,data.FRF.f,data.CO.AsAg),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)'),legend('(Xsi+Xig+Xg)/Xg','As/Ag')


%
% -----Modelo para identificacao n4sid------

%close all

clear SysExp SysIdent Freq Response

global M K fnzId wni wns


%Modelo Experimental
%FRFmod='D';  %FRFs obtidas em termos de Deslocamentos
FRFmod='A';  %FRFs obtidas em termos de Aceleracoes

%cortar o sinal no intervalo pretendido
ffr1=find(data.FRF.f>=frange(1));fr1=ffr1(1);
ffr2=find(data.FRF.f<=frange(2));fr2=ffr2(end);

Freq=data.FRF.f(fr1:fr2)*2*pi;
nFreq=length(Freq);

%FTs: Xig/Ag, Xsg/Ag, Vig/Ag, Vsg/Ag; com: Xsg=Xsi+Xig; Ag=s^2*Xg; V=s*X;
% Response(1,1,:)=data.FRF.XigS2Xg(fr1:fr2);
% Response(2,1,:)=data.FRF.XigS2Xg(fr1:fr2).*complex(zeros(nFreq,1),Freq);
% Response(3,1,:)=data.FRF.XigS2Xg(fr1:fr2)+data.FRF.XsiS2Xg(fr1:fr2);
% Response(4,1,:)=(data.FRF.XigS2Xg(fr1:fr2)+data.FRF.XsiS2Xg(fr1:fr2)).*complex(zeros(nFreq,1),Freq);

if strcmp(FRFmod,'D')
%FTs: Xig/Ag, Xsg/Ag, Ai/Ag, As/Ag; com: Xsg=Xsi+Xig; A=s^2*X; 
Response(1,1,:)=data.FRF.XigS2Xg(fr1:fr2);
Response(2,1,:)=data.FRF.XigS2Xg(fr1:fr2)+data.FRF.XsiS2Xg(fr1:fr2);
Response(3,1,:)=data.FRF.XiXg(fr1:fr2);%data.FRF.XiXg(fr1:fr2);
Response(4,1,:)=data.FRF.XsXg(fr1:fr2);

elseif strcmp(FRFmod,'A')
Response(1,1,:)=data.FRF.XigAg(fr1:fr2)/1000;
Response(2,1,:)=(data.FRF.XigAg(fr1:fr2)+data.FRF.XsiAg(fr1:fr2))/1000;
Response(3,1,:)=data.FRF.AiAg(fr1:fr2);%data.FRF.XiXg(fr1:fr2);
Response(4,1,:)=data.FRF.AsAg(fr1:fr2);
%
end
SysExp=idfrd(Response,Freq,0);

%Modelo Identificado
SysIdent = n4sid(SysExp,1:10,'foc','stability','nk',1,'N4Weight','auto');

[WnId,ZId]=damp(SysIdent.a);fnzId0=sortrows([WnId/2/pi,ZId]);
fnzId=fnzId0([1 3],:);


figure,bode(SysExp(1,1),SysIdent(1,1))
figure,bodeplot(SysExp(2,1),SysIdent(2,1))
figure,bode(SysExp(3,1),SysIdent(3,1))
figure,bode(SysExp(4,1),SysIdent(4,1))


%%
%%%%%%%%%%%%%%%%%%%%%%% Modelo Espacial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(mk,'m')
% 1) identificacao assumindo a massa conhecida como valor certo e invariavel
%massa
mi=1790;
ms=3950;
M=diag([mi ms]);    %matriz de massa


%rigidez: estimativa inicial, valores obtidos em ensaio quase-est�tico
ki0=(17+16.8)*10^3;    %2molas:
ks0=(201.1+205.8)*10^3; %2molas:


%resolucao do problema de optimizacao para acertar as frequencias
x0 = [ki0 ks0];                        % Starting guess
%regiao de projeto
Aop=[];bop=[];Aeqop=[];beqop=[];     %igualdades
lbx=x0/2;ubx=x0*2;             %desigualdades

%
%tolerancia na funcao objetivo
options = optimset('fmincon');
optnew = optimset(options,'Algorithm','interior-point',...
    'TolX',1e-30,'TolCon',1e-10,'TolFun',1e-30,...
    'MaxFunEvals',10^6,'MaxIter',10^4);
[x fvalx,exitflagx] = fmincon(@Stiffness2DOF_fo,x0,Aop,bop,Aeqop,beqop,lbx,ubx,[],optnew);
xconsv=[x./lbx;x./ubx]>=0.99999 & [x./lbx;x./ubx]<=1.00001 %opt nos contrangimentos

%rigidez: resolucao do problema de optimizacao
ks=x(2);            
ki=x(1);            
K=[ki+ks -ks;-ks ks];  %matriz de rigidez

elseif strcmp(mk,'k')
% 2) identificacao assumindo a rigidez conhecida como valor certo e invariavel

%rigidez: valores obtidos em ensaio quase-est�tico
ki=(17+16.8)*10^3;    %2molas:
ks=(201.1+205.8)*10^3; %2molas:
K=[ki+ks -ks;-ks ks];  %matriz de rigidez


%---determinacao das massas---
%massas: estimativa inicial
ms0=4000;    %massa da superstrutura
mi0=1750;     %massa da base

%resolucao do problema de optimizacao para acertar a 1�frequencia
x0 = [mi0 ms0];                        % Starting guess
%regiao de projeto
Aop=[];bop=[];Aeqop=[];beqop=[];     %igualdades
lbx=x0/2;ubx=x0*2;             %desigualdades

%
%tolerancia na funcao objetivo
options = optimset('fmincon');
optnew = optimset(options,'Algorithm','interior-point',...
    'TolX',1e-30,'TolCon',1e-10,'TolFun',1e-30,...
    'MaxFunEvals',10^6,'MaxIter',10^4);
[x fvalx,exitflagx] = fmincon(@Mass2DOF_fo,x0,Aop,bop,Aeqop,beqop,lbx,ubx,[],optnew);
xconsv=[x./lbx;x./ubx]>=0.99999 & [x./lbx;x./ubx]<=1.00001 %opt nos contrangimentos

%massas: resolucao do problema de optimizacao
ms=x(2);            %massa da superstrutura
mi=x(1);            %massa da base
M=diag([mi ms]);    %matriz de massa

end

%frequencias
wni=sqrt(ki/(mi+ms));fni=wni/(2*pi);
wns=sqrt(ks/ms);fns=wns/(2*pi);


%n� pisos equivalente
n=16/fns;

%relacao da massa da base/massa de cada piso
rt=mi/(ms/n);


% ---determinacao do amortecimento-----
%close all

% resolucao do problema de optimizacao para acertar a 1�frequencia

%1) Amortecimento classico
% zetai0=0.05;   %fator de amortecimento do sistema de isolamento 
% zetas0=0.005;   %fator de amortecimento da superstrutura 
% 
% cci=2*(mi+ms)*wni;%amortecimento critico do sist. de isolamento
% ci0=zetai0*cci; %damping
% 
% ccs=2*ms*wns; %amortecimento critico da estrutura
% cs0=zetas0*ccs; %amortecimento da estrutura

% y0 = [ci0 cs0];                        % Starting guess
% lby=y0/2000;uby=y0*20;             %desigualdades
% [y fvaly,exitflagy] = fmincon(@ClassDamp2DOF_fo,y0,Aop,bop,Aeqop,beqop,lby,uby,[],optnew);
% ci=y(1); %damping
% cs=y(2); %amortecimento da estrutura
% C=[ci+cs -cs;-cs cs]; %matriz C

%2) Amortecimento Rayleigh C=alpha*M+beta*K (=Caughey c/2DOF)

[y fvaly,exitflagy] = fmincon(@RayleighDamp2DOF_fo,y0,Aop,bop,Aeqop,beqop,lby,uby,[],optnew);
yconsv=[y./lby;y./uby]>=0.99999 & [y./lby;y./uby]<=1.00001  %opt nos contrangimentos

alpha=y(1);beta=y(2);
C=alpha*M+beta*K;



%vetor da forca aplicada
Gf=[1;0];

%espaco de estados: dx/dt=A*x+B*u+G*ab; y=Ce*x+D*u+H*ab
A=[zeros(2,2) eye(2);-M^-1*K -M^-1*C];G=-[0;0;1;1];B=[0;0;M^-1*Gf];
Ce=[1 0 0 0;0 1 0 0;-M^-1*K -M^-1*C];H=[0;0;0;0];D=[0;0;M^-1*Gf];
eeEp=ss(A,[G B],Ce,[H D]); 
SysEp = idfrd(eeEp,Freq);
[WnEp,ZEp] = damp(A);
fnzEp=sortrows([WnEp([1 3])/2/pi ZEp([1 3])]);


figure,bode(SysExp(1,1),SysIdent(1,1),SysEp(1,1))
figure,bodeplot(SysExp(2,1),SysIdent(2,1),SysEp(2,1))
figure,bode(SysExp(3,1),SysIdent(3,1),SysEp(3,1))
figure,bode(SysExp(4,1),SysIdent(4,1),SysEp(4,1))

error=sqrt((fnzId-fnzEp).^2)./fnzId*100


out=[fnzEp(:,1);fnzEp(:,2)*100;mi;ms;ki;ks;alpha;beta;error(:,1);error(:,2)]; open('out')