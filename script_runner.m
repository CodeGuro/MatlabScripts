linear_bool = false;
useLassoNmat_bool =  false;
n_int = 10;
Alimiter_real = 0.8;
numMatSamples_int = 8;
perturbAmount_real = 0.01;
sigmas_vecreal = [ 0.05 ];
lambdas_vecreal = [];
plotting_bool = false;
mistakeThresh_real = 1E-2;

perturbAmount_vec = [ 0.08:0.01:0.17 0.2 0.3 0.4 0.5 ];
numSampling = 1:72:1000;
numMistakes = nan(14,14);

for i=1:length(numSampling)
    numSamples_int = numSampling(i);
    for j=1:length(perturbAmount_vec)
        perturbAmount_real =  perturbAmount_vec(j);
        numMistakes(i,j) = script_univ( n_int, linear_bool, useLassoNmat_bool, Alimiter_real, numMatSamples_int,...
            perturbAmount_real, numSamples_int, mistakeThresh_real, sigmas_vecreal,...
            lambdas_vecreal, plotting_bool );
    end
end

surf(numSampling, perturbAmount_vec, numMistakes);


