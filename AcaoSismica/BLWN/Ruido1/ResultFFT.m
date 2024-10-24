%visualização da FFT do sinal

load Ruido01a15

time=Ruido01a15(1,:)';
sig=Ruido01a15(2,:)';


[freq sigf]=FFTsig(time,sig);

%FFT: amp & phase
amp=2*abs(sigf);
phase=atan2(imag(sigf),real(sigf));

ti=0;
tf=5;
fi=0;
ff=30;

figure,subplot(2,1,1)
plot(time,sig),xlabel('t (s)'),ylabel('Acel. (g)')
title('Time Signal'),xlim([ti tf])
subplot(2,1,2)
plot(freq,amp),xlabel('f (Hz)'),ylabel('|Acel.| (g)')
title('FFT'),xlim([fi ff])


%Auto-spectra
% AS=2*sqrt(sigf.*conj(sigf));
% 
% 
% figure,subplot(2,1,1)
% plot(time,sig),xlabel('t (s)'),ylabel('Acel. (g)')
% title('Time Signal'),xlim([ti tf])
% subplot(2,1,2)
% plot(freq,AS),xlabel('f (Hz)'),ylabel('|Acel.| (g)')
% title('Auto-Spectra'),xlim([fi ff])
