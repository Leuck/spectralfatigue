% exemplo 2, worked examples do nCode

%% data

load('dados_ex_2.mat');

%% fatigue
vida = spectrallife(psd,sn2,0,3,'dirlik',1)
dano = 1/vida

% RESULTADOS NCODE
% dano = 4.454e-5
% vida = 2.245e4
% rms = 5.891e7
% mean = 0
% irr factor = 0.4362
