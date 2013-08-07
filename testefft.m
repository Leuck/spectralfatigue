%
% Fs = 1000;                    % Sampling frequency
% T = 1/Fs;                     % Sample time
% L = 100000;                     % Length of signal
% t = (0:L-1)*T;                % Time vector
% % Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
% x = 70.*sin(2*pi*50.*t) + 100.*sin(2*pi*120.*t); 
% y = x + 22*randn(size(t));     % Sinusoids plus noise
% figure
% plot(Fs.*t(1:50),y(1:50))
% title('Signal Corrupted with Zero-Mean Random Noise')
% xlabel('time (milliseconds)')
% 
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Y = fft(y,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% PSD_1sided = 2*abs(Y(1:NFFT/2+1));
data = dlmread('PSDansysSyy.csv',';');
PSD_1sided = data(:,2);
f = data(:,1);

matdata = csvread('steelSN.txt');
N = matdata(:,1);
S = matdata(:,2);
Pmat = polyfit(S,log10(N),10);
semilogx(N,S)
xlabel('N')
ylabel('S')

Se = 1.8e7;

% Plot single-sided amplitude spectrum.
figure
plot(f,PSD_1sided) 
title('Single-Sided Amplitude Spectrum of S(t)')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')

% PSD moments
m0 = trapz(f, PSD_1sided);
m1 = trapz(f, f .* PSD_1sided);
m2 = trapz(f, f.^2 .* PSD_1sided);
m4 = trapz(f, f.^4 .* PSD_1sided);

fprintf('\nm0 = %f\nm1 = %f\nm2 = %f\nm4 = %f\n',m0,m1,m2,m4);

RS = linspace(0,Se,500);

% Dirlik's PDF
fm = m1/m0 * sqrt(m2/m4);
gamma = m2/sqrt(m0*m4);
D1 = 2*(fm - gamma^2 )/ (1+gamma^2 );
R = (gamma-fm-D1^2) / (1 - gamma - D1 + D1^2);
D2 = (1 - gamma - D1 + D1^2)/(1 - R);
D3 = 1 - D1 - D2;
Q = 5*(gamma - D1 - D2*R)/ (4*D1);
Z = RS./ (2*sqrt(m0));
DPDF = (D1/Q.*exp(-Z./Q) + D2.*Z./R^2.*exp(-Z.^2./2/R^2) + D3.*Z.*exp(-Z.^2./2) ) ./ (2*sqrt(m0));

% Rayleigh's PDF
RPDF = RS./m0 .* exp(-RS.^2./2./m0);

% Gaussian PDF
GPDF = 1/sqrt(2*pi*m0) * exp(-RS.^2./2/m0);

% Narrow Band PDF
NPDF = RS./(4*m0) .* exp(-RS.^2./8/m0);

figure
plot(RS, DPDF, RS, RPDF, RS, GPDF, RS, NPDF)
title('Different PDF`s for a narrow band signal')
xlabel('\sigma')
ylabel('P')
legend('Dirlik´s PDF','Rayleigh´s PDF','Gaussian PDF')

