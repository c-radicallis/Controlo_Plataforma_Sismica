% Identificacao
%%

clear all 
close all

%ficheiro com FRFs
gain='1';  %04;05;06;07;08;09;1

%--------------------------------------------------------------------
data=load(strcat('2DOF_FreqResp_G',gain));

%gama de frequencias
frange=[0.2 5]; 

%Modelo espacial
%Massa/Rigidez
mk='m';  %'m'-massa conhecida; 'k'-rigidez conhecida

%Amortecimento: 1) alpha & beta ~=0; 2) beta=0;
%y0 = [0.1 0.00001];lby=y0*10^-3;uby=y0*10^3;                 
y0 = [0.1 0];lby=y0*10^-3;uby=y0*10^3; Aeqop=[0 1];beqop=0; %c/beta=0;

%--------------------------------------------------------------------------

%funcao a visualizar
FRFv='all'; %'all'-todos;'AgS2Xg';'XigS2Xg';'XsiS2Xg';'XiXg';'XsXg';'XigAg';'XsiAg';'AiAg';'AsAg';


if strcmp(FRFv,'none')==0
if strcmp(FRFv,'all'),FRFn=fieldnames(data(1).FRF);FRFd=FRFn(2:end);
else FRFd=cellstr(FRFv); end


for j=1:length(FRFd)    
figure('Position',[300 10 600 700]),
for i=1:length(data)
    
    dname=strcat('data(',num2str(i),')');
    hold all
    subplot(311),loglog(eval([dname,'.FRF.f']),abs(eval([dname,'.FRF.',FRFd{j}]))),
    xlim(frange),ylabel('Mag'),title(FRFd{j})
    hold off
    
    hold all
    subplot(312),semilogx(eval([dname,'.FRF.f']),angle(eval([dname,'.FRF.',FRFd{j}]))*180/pi),
    xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
    set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
    hold off
    
    hold all
    subplot(313),semilogx(eval([dname,'.FRF.f']),abs(eval([dname,'.CO.',FRFd{j}]))),
    xlim(frange),ylabel('Coherence'),xlabel('f (Hz)')
    hold off
end
end
end

% -----Modelo para identificacao n4sid------

close all

clear SysExp SysIdent Freq Response

global M K fnzId wni wns


%Modelo Experimental
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


%FTs: Xig/Ag, Xsg/Ag, Ai/Ag, As/Ag; com: Xsg=Xsi+Xig; A=s^2*X; 
Response(1,1,:)=data.FRF.XigS2Xg(fr1:fr2);
Response(2,1,:)=data.FRF.XigS2Xg(fr1:fr2)+data.FRF.XsiS2Xg(fr1:fr2);
% Response(3,1,:)=data.FRF.XiXg(fr1:fr2);%data.FRF.XiXg(fr1:fr2);
% Response(4,1,:)=data.FRF.XsXg(fr1:fr2);
Response(3,1,:)=data.FRF.AiAg(fr1:fr2);%data.FRF.XiXg(fr1:fr2);
Response(4,1,:)=data.FRF.AsAg(fr1:fr2);
%
SysExp=idfrd(Response,Freq,0);

%Modelo Identificado
SysIdent = n4sid(SysExp,1:10,'foc','stability','nk',1,'N4Weight','auto');

[WnId,ZId]=damp(SysIdent.a);fnzId0=sortrows([WnId/2/pi,ZId]);
fnzId=fnzId0([1 3],:);


% figure,bode(SysExp(1,1),SysIdent(1,1))
% figure,bodeplot(SysExp(2,1),SysIdent(2,1))
% figure,bode(SysExp(3,1),SysIdent(3,1))
% figure,bode(SysExp(4,1),SysIdent(4,1))



%%%%%%%%%%%%%%%%%%%%%%% Modelo Espacial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(mk,'m')
% 1) identificacao assumindo a massa conhecida como valor certo e invariavel
%massa
mi=1790;
ms=3950;
M=diag([mi ms]);    %matriz de massa


%rigidez: estimativa inicial, valores obtidos em ensaio quase-estático
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

%rigidez: valores obtidos em ensaio quase-estático
ki=(17+16.8)*10^3;    %2molas:
ks=(201.1+205.8)*10^3; %2molas:
K=[ki+ks -ks;-ks ks];  %matriz de rigidez


%---determinacao das massas---
%massas: estimativa inicial
ms0=4000;    %massa da superstrutura
mi0=1750;     %massa da base

%resolucao do problema de optimizacao para acertar a 1ªfrequencia
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


%nº pisos equivalente
n=16/fns;

%relacao da massa da base/massa de cada piso
rt=mi/(ms/n);


% ---determinacao do amortecimento-----
%close all

% resolucao do problema de optimizacao para acertar a 1ªfrequencia

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

%%
figure,bode(SysExp(1,1),SysEp(1,1))
figure,bodeplot(SysExp(2,1),SysEp(2,1))
figure,bode(SysExp(3,1),SysEp(3,1))
figure,bode(SysExp(4,1),SysEp(4,1))