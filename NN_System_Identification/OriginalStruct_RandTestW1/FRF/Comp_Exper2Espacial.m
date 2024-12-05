%Comparacao com o modelo espacial: dominio da frequencia
%%

clear all 
close all

%ficheiro com FRFs
gain='08';  %05;06;07;08;08_g2
fres='005';    %'001'-0.01Hz;'002'-0.02Hz;'005'-0.05Hz;

%gama de frequencias
frange=[0.15 5]; 

%'..._g2'-com acelerometros +sens
%'_002'-resolucao da FRF em frequencia 002Hz
data=load(strcat('2DOF_FreqResp_G',gain,'_',fres,'hz'));


%------Modelo Experimental-------
%cortar o sinal no intervalo pretendido
ffr1=find(data.FRF.f>=frange(1));fr1=ffr1(1);
ffr2=find(data.FRF.f<=frange(2));fr2=ffr2(end);

Freq=data.FRF.f(fr1:fr2)*2*pi;
nFreq=length(Freq);

%FTs: Xig/Ag, Xsg/Ag, Ai/Ag, As/Ag; com: Xsg=Xsi+Xig; A=s^2*X; 
Response(1,1,:)=data.FRF.XigS2Xg(fr1:fr2);
Response(2,1,:)=data.FRF.XigS2Xg(fr1:fr2)+data.FRF.XsiS2Xg(fr1:fr2);
Response(3,1,:)=data.FRF.XiXg(fr1:fr2);%data.FRF.XiXg(fr1:fr2);
Response(4,1,:)=data.FRF.XsXg(fr1:fr2);
%
SysExp=idfrd(Response,Freq,0);


%-------Modelo espacial------
%massa: media dos valores identificados p/fn com df=0.01Hz
mi=1754;            %massa da base
ms=3427;            %massa da superstrutura
M=diag([mi ms]);    %matriz de massa

%rigidez: valores obtidos em ensaio quase-estático
ki=(17+16.8)*10^3;    %2molas:
ks=(201.1+205.8)*10^3; %2molas:
K=[ki+ks -ks;-ks ks];  %matriz de rigidez

%amortecimento
gf=str2double(gain)/10;
alpha=2.0042*gf^2-3.5643*gf+1.8983;
beta=-0.0009*gf+0.0008;
%
C=alpha*M+beta*K;

%frequencias
wni=sqrt(ki/(mi+ms));fni=wni/(2*pi);
wns=sqrt(ks/ms);fns=wns/(2*pi);

%nº pisos equivalente
n=16/fns;

%relacao da massa da base/massa de cada piso
rt=mi/(ms/n);


%vetor da forca aplicada
Gf=[1;0];

%espaco de estados: dx/dt=A*x+B*u+G*ab; y=Ce*x+D*u+H*ab
A=[zeros(2,2) eye(2);-M^-1*K -M^-1*C];G=-[0;0;1;1];B=[0;0;M^-1*Gf];
Ce=[1 0 0 0;0 1 0 0;-M^-1*K -M^-1*C];H=[0;0;0;0];D=[0;0;M^-1*Gf];
eeEp=ss(A,[G B],Ce,[H D]); 
SysEp = idfrd(eeEp,Freq);

[WnEp,ZEp] = damp(A);
fnzEp=sortrows([WnEp([1 3])/2/pi ZEp([1 3])]);


%---FRFs--------
figure,bode(SysExp(1,1),SysEp(1,1))
figure,bode(SysExp(2,1),SysEp(2,1))
figure,bode(SysExp(3,1),SysEp(3,1))
figure,bode(SysExp(4,1),SysEp(4,1))
