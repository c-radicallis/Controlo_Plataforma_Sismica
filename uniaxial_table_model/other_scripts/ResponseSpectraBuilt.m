function ResponseSpectraBuilt(In,grav)

global m k c

%funcao para construir os espetros de resposta a partir das series de dados
%entrada:
%In-matriz com dados ordenados da seguinte forma:
%   coluna 1 - vetor de tempo
%   coluna 2 a N - vetores de series de aceleracao
%
%grav - gravar dados da integracao
%
%------------------------------------------------------------------
Intime=0;                       %tempo inicial
Sttime=In(end,1);               %tempo final
Indt=In(2,1);                   %discretizacao temporal
[Inlength Inseries]=size(In);   %tamanho da matriz de entrada
%
%--------opcoes de simulacao-------------------------------------------
passo=10^-3;    %passo de integracao
options = simset('FixedStep',passo,'Solver','ode3'); %p. integrcao fixo
%
%-----------frequencias: em bandas de oitava--------------
ke=32; %filtro de oitava (1, 1/3, 1/16)
fni=0.125; %freq inicial
fnf=50;  %freq final
ne=ceil((log(fnf)-log(fni))/log(2)*ke);  %nº de pontos
fn1=[fni fni*2.^((1:1:ne)/ke)]; %vetor de frequencias
% ---------------------------

mdrp=zeros(ne+1,Inseries-1);map=mdrp;
for i=1:ne+1    
    %sistema de simulação LNEC
    m=3750;             %mass
    fn=fn1(i);
    wn=2*pi*fn;         %natural frequency
    zeta=0.05;          %damping ratio (%)
    wd=wn*sqrt(1-zeta^2);
    k=wn^2*m;           %stiffness
    cc=2*m*wn;          %citical damping
    c=zeta*cc;          %damping
    %  
    for j=2:Inseries
        % simulacao do modelo 
        input=[In(:,1),In(:,j)];
        [t, xp, y] = sim('OneDOF_original', [Intime Sttime], options,input);
        %valores maximos
        mdrp(i,j-1)=max(abs(y(:,1)));    %max desl. relativo
        map(i,j-1)=max(abs(y(:,3)));     %max aceleracao
    end
end

%valores medios
mdr=mean(mdrp,2);
mda=mean(map,2);

figure('Position',[100 200 1300 400])
subplot(131),plot(In(:,1),In(:,2:end)),
title('time series'),xlabel('t(s)'),ylabel('a (m/s^2)')
subplot(132),semilogx((1./fn1)',map),xlim([0 4]),grid
xlabel('T (s)'),ylabel('a (m/s^2)'),title('response spectra')
subplot(133),semilogx((1./fn1)',map,(1./fn1)',mda,'-k','LineWidth',2.5),xlim([0 4]),grid
xlabel('T (s)'),ylabel('a (m/s^2)'),title('response spectra')
text(0.5,9,'mean spectra','FontSize',14,'BackgroundColor',[1 1 1])

%save file
if grav==1
    data.T=1./fn1;
    data.a=map;
    data.m=mda;
    data.dr=mdrp;
    data.mdr=mdr;
    
    uisave('data','var1')
end