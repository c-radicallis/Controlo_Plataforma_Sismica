function [f,yf]=FFTsig(t,y)

ta=mean(t(2:end)-t(1:end-1)); %tempo de amostragem
Fs=1/ta; %freq de amostragem

L = length(y);                     % Length of signal
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
yf=Y(1:NFFT/2+1);
end
