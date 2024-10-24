% Determinacao do espetro de resposta para acoes consideradas

clear all, close all

%----Acao de entrada-----------------------------------------------------
%escalonamento da acao
sc=1;
%
%------------10 acoes------------------------
%
%10 acoes EC8 DNA Zona 1.1 Solo D Categ. II
% load EC8_A1.txt
% In=EC8_A1;In(:,2:end)=EC8_A1(:,2:end)*sc;
%
%10 acoes EC8 DNA Zona 2.1 Solo D Categ. II
load EC8_A2.txt
In=EC8_A2;In(:,2:end)=EC8_A2(:,2:end)*sc;
%
%-----------1 acao--------------------------------
%
%acao EC8 DNA Zona 1.1 Solo D Categ. II
% load EC8DNA11DIIacel.txt  %sinal sismo - valores em m/s2
% In(:,1)=EC8DNA11DIIacel(:,1);In(:,2)=EC8DNA11DIIacel(:,2)*sc;  
%
%acao EC8 DNA Zona 2.1 Solo D Categ. II
% load EC8DNA21DIIacel.txt  %sinal sismo - valores em m/s2
% In(:,1)=EC8DNA21DIIacel(:,1);In(:,2)=EC8DNA21DIIacel(:,2)*sc;
%
%acao EC8
% load SinalAcel.txt  %sinal sismo - valores em g
% In(:,1)=SinalAcel(:,1);In(:,2)=SinalAcel(:,2)*9.86*sc;
%
%------------------------------------------------------------------
Intime=0;                       %tempo inicial
Sttime=In(end,1);               %tempo final
Indt=In(end,1)-In(end-1,1);     %discretizacao temporal
[Inlength Inseries]=size(In);   %tamanho da matriz de entrada
%
%--------opcoes de simulacao-------------------------------------------
passo=10^-3;    %passo de integracao
options = simset('FixedStep',passo,'Solver','ode3'); %p. integrcao fixo
options2 = simset('Solver','ode45');                %p. integrcao variavel
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

%% comparacao com o espetro pretendido

ac=2;

%-----------espetro de resposta-----------------------------------------
%-----parametros---------------------------------------------------------
ag=2.5;     %aceleracao a superficie
S=1.5;      %coeficiente do solo
na=1;       %coef. correc. amortecimento - na=1 p/zeta=5%
T0=0;   %valor do periodo de inicio
Tb=0.1;     %limite inf. do patamar cte
Tc1=0.8;    %limite sup. do patamar cte (acao 1)
Tc2=0.3;    %limite sup. do patamar cte (acao 2)
Tc=(ac==1)*Tc1 + (ac==2)*Tc2;   %limite sup. do patamar cte
Td=2;       %valor que define o ramo de desl. cte
Tf=8;       %valor do periodo final

%espetro de resposta pretendido
Tcd=(Tc:0.0001:Td)';Tdf=(Td:0.0001:Tf)';
Tp=[T0;Tb;Tc;Tcd;Tdf];
Se=[S*ag;2.5*S*na*ag;2.5*S*na*ag;ag*S*na*2.5*Tc./Tcd;ag*S*na*2.5*Tc*Td./Tdf.^2];
%-------------------------------------------------------------------------

% figure('Position',[100 200 1300 400])
% subplot(131),plot(tempo,acel),
% title('time series'),xlabel('t(s)'),ylabel('a (m/s^2)')
% subplot(132),plot(Tp,Se,(1./fn1)',map),xlim([0 4]),grid
% xlabel('T (s)'),ylabel('a (m/s^2)'),title('response spectra')
% subplot(133),plot((1./fn1)',mda,Tp,Se,'-r','LineWidth',1.5),xlim([0 4]),grid
% xlabel('T (s)'),ylabel('a (m/s^2)'),title('response spectra')
% legend('mean spectra','EC8 spectra')

%log scale
figure('Position',[100 200 500 400])
semilogx((1./fn1)',mda,Tp,Se,'-r','LineWidth',1.5),xlim([0 4]),grid
xlabel('T (s)'),ylabel('a (m/s^2)'),title('response spectra')
legend('mean spectra','EC8 spectra')

%linear scale
figure('Position',[610 200 500 400])
plot((1./fn1)',mda,Tp,Se,'-r','LineWidth',1.5),xlim([0 4]),grid
xlabel('T (s)'),ylabel('a (m/s^2)'),title('response spectra')
legend('mean spectra','EC8 spectra')

%% visualiuzar espetro especifico
close all

for i=1:Inseries
espetro=i;

figure,plot((1./fn1)',map(:,espetro)),xlim([0 4]),grid
xlabel('T (s)'),ylabel('a (m/s^2)'),title('response spectra')
end
