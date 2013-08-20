% exemplo 2, worked examples do nCode

%% data

load('dados_ex_2.mat');

%% fatigue
vida(1) = spectrallife(psd,sn1,0,3,'Dirlik',0);
vida(2) = spectrallife(psd,sn2,0,3,'Dirlik',1);
vida(3) = spectrallife(psd,sn3,0,3,'Dirlik',0);

vidamedia = mean(vida)
dano = 1/vidamedia

% RESULTADOS NCODE
% dano = 4.454e-5
% vida = 2.245e4
% rms = 5.891e7
% mean = 0
% irr factor = 0.4362
