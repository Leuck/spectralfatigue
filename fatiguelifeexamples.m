% Vida em fadiga, diversos exemplos

%% Dados de entrada
sinal = load('sinal_baixa_amp.mat');
sinal2 = load('sinal_campo.mat');
t=linspace(0,20,20*1300)';
teste1.t = t;
teste2.t = t;
teste1.ext=277.92/2*sin(2*pi*5*t);
teste2.ext=253.216/2*sin(2*pi*5*t);

sn = csvread('1006SN.txt');
pdf = 'dirlik'; % usar dirlik, rayleigh, gauss ou narrow
showplots = 0; % =1 para mostrar gráficos

%% Chamadas de funcao
for criteria = [1 2 3]
    T(1) = spectralfatigue(teste1,sn,criteria,pdf,0); % 5Hz, 277.92MPa -> 111h
    T(2) = spectralfatigue(teste2,sn,criteria,pdf,0); % 5Hz, 253.216MPa -> 333h
    T(3) = spectralfatigue(sinal2,sn,criteria,pdf,1); % rainflow --> 1152h
    %T(4) = spectralfatigue(sinal,sn,criteria,pdf,showplots); % Tf analitico = 735h
        
    %% Vida em fadiga
    fprintf('Criterio %d:\n',criteria)
    Tv = T/3600;
    for i = Tv
        fprintf('%.0f\t',i)
    end
    fprintf('\n\n')
end