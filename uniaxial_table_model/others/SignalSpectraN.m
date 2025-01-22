function SignalSpectraN(sinal,le,ov,tw,tp,spec)

%funcao para ver o sinal no tempo com valor RMS e o 
%espectro de Welch do sinal
%entradas: 
%>matriz com o sinal: coluna 1-vetor de tempos; 
%                     colunas 2 a N - respostas
%
%parametros do estimador espectral
%le: comprimento do vector de tempo para determinar o espectro de Welch
%ov: perentagem de sobreposicao entre segmentos
%tw: tukey window type: 0 to 1; 0={rectangular};1={hanning}
%tp: escala do espectro: 1-linear; 2-logaritmica
%spec: espetro em MSS (Mean Squared Spectrum) ou PSD (Power Spectral Sensity)
%
%-------------------------------------------------

[slength sseries]=size(sinal);   %tamanho da matriz de entrada
tlim=sinal(end,1);               %tempo final

dt=sinal(end,1)-sinal(end-1,1); %disc no tempo
Fs=1/dt;        %freq amostragem

%Welch's averaged, modified periodogram spectral estimation
hp=spectrum.welch({'tukey',tw},le,ov);  %{'tukey',1}={hanning}
%hp=spectrum.periodogram;

for i=2:sseries
clear ms

if spec==1
%Mean-square (power) spectrum
ms = msspectrum(hp,sinal(:,i),'Fs',Fs); 
%Freq e amplitude
msd(:,i-1)=2*sqrt(ms.Data/2);

else
ms=psd(hp,sinal(:,i),'Fs',Fs);
msd(:,i-1)=ms.Data;
end

end

msf=ms.Frequencies;

amax=max(max(sinal(:,2:end)));
amin=min(min(sinal(:,2:end)));

figure,subplot(2,1,1)
plot(sinal(:,1),sinal(:,2:end),[0 tlim],[amax amax],'--k',...
    [0 tlim],[amin amin],'--k'),
xlim([0 tlim]),xlabel('t (s)'),ylabel('A (unit)')
title('time signal'),
legend('signal',['min value =',sprintf('%5.4f',amin)],...
    ['Max value =',sprintf('%5.4f',amax)])
subplot(2,1,2)
if tp==1
    plot(msf,msd),
elseif tp==2
    loglog(msf,msd),
end
xlabel('f (Hz)')
if spec==1
    ylabel('|A| (unit)')
    title('Fourier Transform')
else
    ylabel('|A| (unit^2/Hz)')
    title('Power Spectral Density')
end
