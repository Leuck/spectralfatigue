% One sided PSD of a time signal
function [PSD_1sided, f] = psd1s(sinal)
Fs = inv(max(sinal.t)/length(sinal.t));	% Sampling frequency
L = length(sinal.t);                    % Length of signal
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(sinal.ext,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1)';
PSD_1sided = 2*abs(Y(1:NFFT/2+1));

% size(PSD_1sided)
% size(f)

% Fs = round(inv(max(sinal.t)/length(sinal.t)));
% [PSD_1sided,f] = pwelch(sinal.ext,[],[],[],Fs,'onesided');
% PSD_1sided = sqrt(PSD_1sided).*(f+.0000000001);