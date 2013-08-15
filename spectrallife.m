% Fatigue life estimate - spectral method
% sn should be of the form [N S]
% Ricardo Frederico Leuck Filho 2012/2-2013/1
function Tf = spectrallife(psd,sn,meanstress,criteria,pdf,showplots)
% %% Options
% criteria = 2; % Mean stress correction (1=Goodman, 2=Gerber, 3=sem correcao)
% showplots = 0; % Show plots? 0=no, 1=yes

f = psd(:,1);
PSD_1sided = psd(:,2);

%% PSD moments
frad=f;%*2*pi;
% size(frad)
m0 = trapz(frad, PSD_1sided);
m1 = trapz(frad, frad .* PSD_1sided);
m2 = trapz(frad, frad.^2 .* PSD_1sided);
m4 = trapz(frad, frad.^4 .* PSD_1sided);
fprintf('\nMoments:\t%.4E\t%.4E\t%.4E\t%.4E\n',m0,m1,m2,m4);
umsigma = sqrt(m0)

%% S-N curve with Mean Stress Correction
N = sn(:,1);
S = sn(:,2); clear sn;

Sm = meanstress;	% mean stress
Su = max(S);
if criteria == 1
    Sc = S * (1-Sm/Su);        % Goodman
    fprintf('\nMean Stress Correction: Goodman')
elseif criteria == 2    
    Sc = S * (1-(Sm/Su)^2);    % Gerber
    fprintf('\nMean Stress Correction: Gerber')
else
    Sc = S;
end
%% Stress Range
% No extrapolation, the stress range is given by the S-N curve
Smax = max(Sc);
Smin = min(Sc); % fatigue limit
CRS = linspace(Smin,Smax,500); % range for calculation

%% Fatigue life estimation
% Expected damage rate
NF=Ncycles(N,Sc, CRS);
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
    % Plot single-sided amplitude spectrum.
    subplot(2,1,1)
    semilogy(f,PSD_1sided) 
    title('Single-Sided Amplitude Spectrum of S(t)'); xlabel('Frequency, Hz')
    ylabel('S(f)^2/f')
    
    % Plot S-N Diagram
    subplot(2,2,3)
    semilogx(N,S,'*-',N,Sc,N,Sc); xlabel('N, cycles'); ylabel('S, Pa')
    title('S-N Curve'); 

    % Plot Different PDF's
    subplot(2,2,4)
    plot(CRS, PDF)
    title(strcat({pdf},' PDF')); xlabel('\sigma, Pa'); ylabel('P')
end

%% Expected Life
Tf = 1/Ed;
if showplots ~= 0
    fprintf('\nVida: %.0f horas.\n',Tf/3600)
end
