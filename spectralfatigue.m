% Fatigue life estimate - spectral method
% time vector must be in sinal.t, signal must be in sinal.ext
% sn should be of the form [N S]
% Ricardo Frederico Leuck Filho 2012/2-2013/1
function Tf = spectralfatigue(sinal,sn,criteria,pdf,showplots)
% %% Options
% criteria = 2; % Mean stress correction (1=Goodman, 2=Gerber, 3=sem correção)
% showplots = 0; % Show plots? 0=no, 1=yes

%% Stress PSD
Fs = inv(max(sinal.t)/length(sinal.t));	% Sampling frequency 1200; %
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

%% PSD moments
frad=f*2*pi;
% size(frad)
m0 = trapz(frad, PSD_1sided);
m1 = trapz(frad, frad .* PSD_1sided);
m2 = trapz(frad, frad.^2 .* PSD_1sided);
m4 = trapz(frad, frad.^4 .* PSD_1sided);
%fprintf('\nm0 = %.4E\nm1 = %.4E\nm2 = %.4E\nm4 = %.4E',m0,m1,m2,m4);

%% S-N curve with Mean Stress Correction
N = sn(:,1);
S = sn(:,2); clear sn;
%S = S/1e6;
%S = S*6.895;

Sm = mean(sinal.ext);	% mean stress
Su = 500;               % Ultimate tensile strength of SAE 1008
Sc(:,3) = S;
Sc(:,1) = S * (1-Sm/Su);        % Goodman
Sc(:,2) = S * (1-(Sm/Su)^2);    % Gerber

%% Stress Range
% No extrapolation, the stress range is given by the S-N curve
Smax = max(Sc(:,criteria));
Smin = min(Sc(:,criteria)); % fatigue limit
CRS = linspace(Smin,Smax,500); % range for calculation

%% Fatigue life estimation
% Expected damage rate
NF=Ncycles(N,Sc(:,criteria), CRS);
switch lower(pdf)
    case {'dirlik'}
        % Dirlik's PDF
        xm = m1/m0 * sqrt(m2/m4);
        mu = sqrt(m4/m2);
        sigmax = sqrt(m0);
        lambda0 = sqrt(m2/m0);
        gamma = lambda0/mu;
        z = abs(CRS./ (2*sigmax));
        C1 = 2*(xm - gamma^2 )/ (1+gamma^2 );
        alpha = (gamma-xm-C1^2) / (1 - gamma - C1 + C1^2);
        C2 = (1 - gamma - C1 + C1^2)/(1 - alpha);
        C3 = 1 - C1 - C2;
        tau = 1.25*(gamma - C3 - C2*alpha)/ C1;
        PDF = (C1/tau.*exp(-CRS/2/sigmax./tau) + C2.*(CRS/2/sigmax)./alpha^2.*exp(-(CRS/2/sigmax).^2./2/alpha^2) + C3.*(CRS/2/sigmax).*exp(-(CRS/2/sigmax).^2./2) );

        Ed = mu.*trapz( CRS, PDF./NF);
    case {'rayleigh'}
        % Rayleigh's PDF
        PDF = CRS./m0 .* exp(-CRS.^2./2./m0);
        Ed = sqrt(m4/m2).*trapz( CRS, PDF./NF);
    case {'gauss'}
        % Gaussian PDF
        PDF = 1/sqrt(2*pi*m0) * exp(-CRS.^2./2/m0);
        Ed = sqrt(m4/m2).*trapz( CRS, PDF./NF);
    case {'narrow'}
        % Narrow Band PDF
        PDF = CRS./(4*m0) .* exp(-CRS.^2./8/m0);
        Ed = sqrt(m4/m2).*trapz( CRS, PDF./NF);
    otherwise
        fprintf('\nError: pdf should be a string. Either dirlik, rayleigh, gauss or narrow.\n')
end

%% Plot Graphs
if showplots ~= 0
    % Plot time history
    figure;
    subplot(2,2,1)
    plot(sinal.t(1:300),sinal.ext(1:300));
    title('Time History Sample'); xlabel('Time, s'); ylabel('Stress, MPa')
    
    % Plot single-sided amplitude spectrum.
    subplot(2,2,2)
    plot(f,PSD_1sided) 
    title('Single-Sided Amplitude Spectrum of S(t)'); xlabel('Frequency, Hz')
    ylabel('|S(f)|, MPa')
    
    % Plot S-N Diagram
    subplot(2,2,3)
    semilogx(N,S,'*-',N,Sc(:,1),N,Sc(:,2)); xlabel('N, cycles'); ylabel('S, MPa')
    title('S-N Curve'); legend('Zero Mean','Goodman','Gerber')

    % Plot Different PDF's
    subplot(2,2,4)
    plot(CRS, PDF)
    title(strcat({pdf},' PDF')); xlabel('\sigma, MPa'); ylabel('P')
end

%% Expected Life
Tf = 1/Ed;
if showplots ~= 0
    fprintf('\nVida: %.0f horas.\n',Tf/3600)
end
