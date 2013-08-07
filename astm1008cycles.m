% Computes number of cycles to failure for a given stress.
function Ni = astm1008cycles(Si)

% Si is the stress amplitude for which the life in cycles will be computed. 
B = -0.073;
A = 513.36;
Ni = (Si/A).^(1/B);
