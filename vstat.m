function vstat(A)
M = mean(A);
S = std(A);
SI = size(A);
MA = max(A);
mA = min(A);

fprintf('\nMean\t\tStd\t\t\tfrom\t\tto\t\t\tsize\n%.2E\t%.2E\t%.2E\t%.2E',M,S,mA,MA)
fprintf('\t%d by %d\n\n',SI(1),SI(2))