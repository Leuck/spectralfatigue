
sinal = load("sinal_campo.mat");

Fs = inv(max(sinal.t)/length(sinal.t));  % Sampling frequency
T = 1/Fs;                     % Sample time
L = length(sinal.t);          % Length of signal
t = (0:L-1)*T;                % Time vector

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(sinal.ext,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1)';
PSD_1sided = 2*abs(Y(1:NFFT/2+1));

figure
plot(f,PSD_1sided);
title('PSD Unilateral das tensões atuantes')
xlabel('f(Hz)')

