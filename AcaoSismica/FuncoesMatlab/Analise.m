clear
close all

%sinal
Fs=1000;
t=0:1/Fs:.3;
nt=length(t);
x=cos(2*pi*t*200)+randn(size(t));

%espectro
Hs=spectrum.periodogram;	     % Periodogram spectrum
%Hs = spectrum.welch; %Welch's averaged, modified periodogram

% mean-squared spectrum (MSS)
% The MSS of a signal is the Fourier transform of that signal's
% autocorrelation.
msspect = msspectrum(Hs,x,'Fs',Fs);

%transformada de fourier
[f,yf]=FFTsig(t,x);
yf1=yf.*conj(yf);  %valor x conjugado
yf1m=2*sqrt(yf1);  %one-sided
yf2=2*abs(yf);
yf3=2*abs(yf).^2;
yf4=2*yf1;         %one-sided auto-spectra     

%comparacao
figure,semilogy(f,(yf1m),f,(yf2),f,(yf3),...
    f,(yf4))
legend('2 sqrt(yf yf*)','2 |yf|','2 |yf|^2','2 yf yf*')


%outra possibilidade de compor o MSS
Hmss=dspdata.msspectrum(2*yf1,'Fs',Fs,'spectrumtype','onesided'); 


%comparacao de espectros
figure,semilogy(msspect.Frequencies,(msspect.Data),f,(yf4),...
    Hmss.Frequencies,(Hmss.Data))
legend('msspec1','via fft','mssspec2')
xlabel('f (Hz)'),ylabel('Power (amp^2)')

figure,semilogy(msspect.Frequencies,2*sqrt(msspect.Data/2),f,(yf1m))
title('Fourier Transform'),legend('via msspec','via fft')

%Power Spectral density
psddata=psd(Hs,x,'Fs',Fs);

nfft = 2^nextpow2(nt);
Pxx = abs(fft(x,nfft)).^2/length(t)/Fs;
freq = Fs/2*linspace(0,1,nfft/2+1);
Pxxf=2*Pxx(1:nfft/2+1);

% Create a single-sided spectrum
Hpsd = dspdata.psd(Pxxf,'Fs',Fs);  
%figure,plot(Hpsd); 

%outra forma de construir a PSD 
[Pxx2,f2]=periodogram(x,[],'onesided',nfft,Fs);
[Pxx2n,f2n]=periodogram(x,[],'onesided',nfft);
f3=f2n*(Fs/2)/pi; Pxx3=Pxx2n*pi/(Fs/2); %normalizado
figure,semilogy(f2,Pxx2,f3,Pxx3),title('Power Spectral Density') 
legend('via non-norm','via norm')
xlabel('f (Hz)'),ylabel('Power/frequency (amp^2/Hz)')

%comparacao PSD
figure,semilogy(psddata.Frequencies,(psddata.Data),...
    Hpsd.Frequencies,(Hpsd.Data),freq,(Pxxf),f2,Pxx2)
title('Power Spectral Density')
legend('via psd','via dsp.data','via fft','via periodogram')
xlabel('f (Hz)'),ylabel('Power/frequency (amp^2/Hz)')


%comparacao com o calculo via fft
figure,plot(freq,10*log10(Pxxf),f,10*log10(yf4*nt/Fs))
title('Power Spectral Density')
legend('via fft','via fft/L')
xlabel('f (Hz)'),ylabel('Power/frequency (amp^2/Hz)')

%figure,semilogy(freq,(Pxxf),msspect.Frequencies,(msspect.Data)*nt/Fs)
