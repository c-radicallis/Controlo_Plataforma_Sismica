function Integra(In,corr,grav)

%funcao para visualizar a integracao do sinal de entrada, o qual esta
%representado em termos de aceleracoes

%entrada:
%In - matriz data com dados ordenados da seguinte forma:
%   coluna 1 - vetor de tempo
%   coluna 2 a N - vetores de series de aceleracao
%
%corr - parametro para correcao do deslocamento
%
%grav - gravar dados da integracao
%
%------------------------------------------------------------------------
%
Intime=0;                       %tempo inicial
Sttime=In(end,1);               %tempo final
Indt=In(end,1)-In(end-1,1);     %discretizacao temporal
[Inlength Inseries]=size(In);   %tamanho da matriz de entrada

velc=zeros(Inlength,Inseries-1);xc=velc;acelc=xc;xc1=xc;
for i=2:Inseries
    acelc(:,i-1)=In(:,i);    %entrada: aceleracao
    for j=1:Inlength-1
        velc(j+1,i-1)=velc(j,i-1)+acelc(j,i-1)*(In(j+1,1)-In(j,1)); %velocidade
        xc1(j+1,i-1)=xc1(j,i-1)+velc(j,i-1)*(In(j+1,1)-In(j,1));      %deslocamento
    end
    if corr==1
        k=1;
        while abs(xc1(k,i-1))<10^-3
            k=k+1;
        end
        m=(xc1(end,i-1)-xc1(k,i-1))/(Sttime-In(k,1));
        xc(1:k-1,i-1)=xc1(1:k-1,i-1);
        xc(k:end,i-1)=xc1(k:end,i-1)-m*(In(k:end,1)-In(k,1));
    else
        xc(:,i-1)=xc1(:,i-1);
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

%save file
if grav==1
    data.time=In(:,1);
    data.x=xc;
    data.v=velc;
    data.a=acelc;
    uisave('data','var1')
end