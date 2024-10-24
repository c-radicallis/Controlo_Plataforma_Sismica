%% Visualizacao do sinal no tempo - comp c/simulink

clear all, %close all

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
% load EC8_A2.txt
% In=EC8_A2;In(:,2:end)=EC8_A2(:,2:end)*sc;
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
load SinalAcel.txt  %sinal sismo - valores em g
In(:,1)=SinalAcel(:,1);In(:,2)=SinalAcel(:,2)*9.86*sc;
%
%------------------------------------------------------------------
Intime=0;                       %tempo inicial
Sttime=In(end,1);               %tempo final
Indt=5*10^-3;%In(end,1)-In(end-1,1);     %discretizacao temporal
[Inlength Inseries]=size(In);   %tamanho da matriz de entrada
%
%--------opcoes de simulacao-------------------------------------------
passo=Indt;    %passo de integracao
options = simset('FixedStep',passo,'Solver','FixedStepDiscrete'); %p. integrcao fixo
%
%-------------------------------------------------------------------
%% Integracao


velc=zeros(Inlength,Inseries-1);xc=velc;acelc=xc;
for i=2:Inseries
    acelc(:,i-1)=In(:,i);    %entrada: aceleracao
    for j=1:Inlength-1
        velc(j+1,i-1)=velc(j,i-1)+acelc(j,i-1)*(In(j+1,1)-In(j,1)); %velocidade
        xc(j+1,i-1)=xc(j,i-1)+velc(j,i-1)*(In(j+1,1)-In(j,1));      %deslocamento
    end
end

%valores máximos e minimos
xmax=max(max(xc));
xmin=min(min(xc));

vmax=max(max(velc));
vmin=min(min(velc));

amax=max(max(acelc));
amin=min(min(acelc));

%
%apresentacao de resultados
figure,
subplot(311),plot(In(:,1),acelc,[0 Sttime],[amax amax],'--k',...
    [0 Sttime],[amin amin],'--k'),ylabel('a (m/s^2)'),title('time series')
text(Sttime,amax,num2str(amax,'%5.4f'),'HorizontalAlignment','left')
text(Sttime,amin,num2str(amin,'%5.4f'),'HorizontalAlignment','left')
subplot(312),plot(In(:,1),velc,[0 Sttime],[vmax vmax],'--k',...
    [0 Sttime],[vmin vmin],'--k'),ylabel('v (m/s)')
text(Sttime,vmax,num2str(vmax,'%5.4f'),'HorizontalAlignment','left')
text(Sttime,vmin,num2str(vmin,'%5.4f'),'HorizontalAlignment','left')
subplot(313),plot(In(:,1),xc,[0 Sttime],[xmax xmax],'--k',...
    [0 Sttime],[xmin xmin],'--k'),ylabel('x (m)'),xlabel('t (s)')
text(Sttime,xmax,num2str(xmax,'%5.4f'),'HorizontalAlignment','left')
text(Sttime,xmin,num2str(xmin,'%5.4f'),'HorizontalAlignment','left')


%% Filtragem

%filtrar sinais
% 
Xout=zeros(size(In));Xout(:,1)=In(:,1);
for i=2:Inseries
    % filtragem
    input=[In(:,1),In(:,i)];
    %[t, xp, y] = sim('Filter015', [Intime Sttime], options,input);
    [t, xp, y] = sim('FilterLowPass50Hz', [Intime Sttime], options,input);
    Xout(:,i)=y;
end


velc=zeros(Inlength,Inseries-1);xc=velc;acelc=xc;
for i=2:Inseries
    acelc(:,i-1)=Xout(:,i);    %entrada: aceleracao
    for j=1:Inlength-1
        velc(j+1,i-1)=velc(j,i-1)+acelc(j,i-1)*(In(j+1,1)-In(j,1)); %velocidade
        xc(j+1,i-1)=xc(j,i-1)+velc(j,i-1)*(In(j+1,1)-In(j,1));      %deslocamento
    end
end


%valores máximos e minimos
xmax=max(max(xc));
xmin=min(min(xc));

vmax=max(max(velc));
vmin=min(min(velc));

amax=max(max(acelc));
amin=min(min(acelc));

%
%apresentacao de resultados
figure,
subplot(311),plot(In(:,1),acelc,[0 Sttime],[amax amax],'--k',...
    [0 Sttime],[amin amin],'--k'),ylabel('a (m/s^2)'),title('time series')
text(Sttime,amax,num2str(amax,'%5.4f'),'HorizontalAlignment','left')
text(Sttime,amin,num2str(amin,'%5.4f'),'HorizontalAlignment','left')
subplot(312),plot(In(:,1),velc,[0 Sttime],[vmax vmax],'--k',...
    [0 Sttime],[vmin vmin],'--k'),ylabel('v (m/s)')
text(Sttime,vmax,num2str(vmax,'%5.4f'),'HorizontalAlignment','left')
text(Sttime,vmin,num2str(vmin,'%5.4f'),'HorizontalAlignment','left')
subplot(313),plot(In(:,1),xc,[0 Sttime],[xmax xmax],'--k',...
    [0 Sttime],[xmin xmin],'--k'),ylabel('x (m)'),xlabel('t (s)')
text(Sttime,xmax,num2str(xmax,'%5.4f'),'HorizontalAlignment','left')
text(Sttime,xmin,num2str(xmin,'%5.4f'),'HorizontalAlignment','left')



%% Comparacao com e sem filtro

%serie nº
se=4;


%
velo=zeros(Inlength,1);xo=velo;
acelo=In(:,se+1);    %entrada: aceleracao
for j=1:Inlength-1
     velo(j+1)=velo(j)+acelo(j)*(In(j+1)-In(j)); %velocidade
     xo(j+1)=xo(j)+velo(j)*(In(j+1)-In(j));      %deslocamento
end


Xout=zeros(Inlength,2);Xout(:,1)=In(:,1);
% filtragem
input=[In(:,1),In(:,se+1)];
%[t, xp, y] = sim('Filter015', [Intime Sttime], options,input);
[t, xp, y] = sim('FilterLowPass50Hz', [Intime Sttime], options,input);
Xout(:,2)=y;



velc=zeros(Inlength,1);xc=velc;
acelc=Xout(:,2);    %entrada: aceleracao
for j=1:Inlength-1
     velc(j+1)=velc(j)+acelc(j)*(In(j+1)-In(j)); %velocidade
     xc(j+1)=xc(j)+velc(j)*(In(j+1)-In(j));      %deslocamento
end

%
%apresentacao de resultados
figure,
subplot(311),plot(In(:,1),acelo,In(:,1),acelc),
ylabel('a (m/s^2)'),title('time series')
subplot(312),plot(In(:,1),velo,In(:,1),velc),ylabel('v (m/s)')
subplot(313),plot(In(:,1),xo,In(:,1),xc),ylabel('x (m)'),xlabel('t (s)')
legend('original','filtradao')

%% Visualizacao do sinal no tempo - comp c/simulink

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
passo=5*10^-3;    %passo de integracao
options = simset('FixedStep',passo,'Solver','ode3'); %p. integrcao fixo
options2 = simset('Solver','ode45');                %p. integrcao variavel
%
%--------------------------------------------------------------------------
%
velc=zeros(Inlength,Inseries-1);xc=velc;acelc=xc;
Xin=zeros(size(In));Vin=Xin;
for i=2:Inseries
    acelc(:,i-1)=In(:,i);    %entrada
    % simulacao do modelo
    input=[In(:,1),In(:,i)];
    [t, xp, y] = sim('IntegA', [Intime Sttime], options,input);
    Vin(:,i)=y(:,1);
    Xin(:,i)=y(:,2);
    for j=1:Inlength-1
        velc(j+1,i-1)=velc(j,i-1)+acelc(j,i-1)*(In(j+1,1)-In(j,1));
        xc(j+1,i-1)=xc(j,i-1)+velc(j,i-1)*(In(j+1,1)-In(j,1));
    end
end


figure,
subplot(311),plot(In(:,1),In(:,2:end)),ylabel('a (m/s^2)')
subplot(312),plot(In(:,1),velc,t,Vin(:,2:end)),ylabel('v (m/s)')
subplot(313),plot(In(:,1),xc,t,Xin(:,2:end)),ylabel('x (m)'),xlabel('t (s)')

