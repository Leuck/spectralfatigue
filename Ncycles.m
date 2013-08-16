% Computes number of cycles to failure for a given stress.
function Ni = Ncycles(N,S,Si)
% N and S are the material's S-N diagram.
% Si is the stress amplitude for which the life in cycles will be computed. 
LS = log10(S);
LN = log10(N);
LSi = log10(Si);

P = polyfit(LS,LN,2);
LNi = polyval(P,LSi);
Ni = 10.^LNi;
