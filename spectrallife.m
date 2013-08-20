% Fatigue life estimate - spectral method
% psd should be of the form [frequency spectrum]
% sn should be of the form [N S]
% Ricardo Frederico Leuck Filho 2012/2-2013/1
function Tf = spectrallife(psd,sn,meanstress,criteria,pdf,showplots)
% %% Options
% meanstress: mean stress
% criteria: Mean stress correction (1=Goodman, 2=Gerber, 3=no mean value correction)
% pdf: Probability density function estimator, 'dirlik', 'gauss',
% 'rayleigh' or 'narrow'.
% showplots: Show plots? 0=no, 1=yes

nop = 200;      % number of points for calculations
f = psd(:,1);   % frequency
PSD_1sided = psd(:,2);  % response stress psd at point of interest

%% PSD moments
m0 = trapz(f, PSD_1sided);
m1 = trapz(f, f .* PSD_1sided);
m2 = trapz(f, f.^2 .* PSD_1sided);
m4 = trapz(f, f.^4 .* PSD_1sided);
if showplots ~= 0
    %fprintf('\nMoments:\t%.4E\t%.4E\t%.4E\t%.4E\n',m0,m1,m2,m4);
    fprintf('\nRMS: %.4E\tIrr: %.4E\tPeaks: %.4E\n',sqrt(m0),sqrt((m2^2)/(m0*m4)),sqrt(m4/m2));
end
%% S-N curve with Mean Stress Correction
N = sn(:,1);
S = sn(:,2);

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
CRS = linspace(Smin,Smax,nop); % range for calculation

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
        z = abs(CRS./(2*sigmax));
        C1 = 2*(xm - gamma^2 )/ (1+gamma^2 );
        alpha = (gamma-xm-C1^2) / (1 - gamma - C1 + C1^2);
        C2 = (1 - gamma - C1 + C1^2)/(1 - alpha);
        C3 = 1 - C1 - C2;
        tau = 1.25*(gamma - C3 - C2*alpha)/ C1;
        PDF1 = C1/tau.*exp(-z/tau);
        PDF2 = (C2*z/alpha^2).*exp(-(z.^2)./(2*alpha^2));
        PDF3 = C3.*z.*exp(-(z.^2)./2);
        PDF = (PDF1 + PDF2 + PDF3)./(2*sigmax);
        integrand = PDF./NF;
        Ed = mu.*trapz( CRS, integrand);
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
    subplot(2,2,1)
    semilogy(f,PSD_1sided) 
    title('Single-Sided Amplitude Spectrum of S(t)'); xlabel('Frequency, Hz')
    ylabel('S(f)^2/f')
    
    % Plot S-N Diagram
    subplot(2,2,2)
    semilogx(N,S,'*-',N,Sc,N,Sc); xlabel('N, cycles'); ylabel('S, Pa')
    title('S-N Curve'); 

    % Plot Different PDF's
    subplot(2,2,3)
    semilogy(CRS, NF)
    title('N'); xlabel('\sigma, Pa'); ylabel('N')

    % Plot Different PDF's
    subplot(2,2,4)
    plot(CRS, PDF)
    title(strcat({pdf},' PDF')); xlabel('\sigma, Pa'); ylabel('P')
end

%% Expected Life
Tf = 1/Ed;
if showplots ~= 0
    fprintf('\nVida: %.2f horas. Taxa de dano: %.3E\n',Tf/3600,Ed)
end
