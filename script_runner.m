linear_bool = false;
useLassoNmat_bool =  false;
n_int = 10;
Alimiter_real = 0.9;
numMatSamples_int = 8;
sigmas_vecreal = [ 0.05 ];
lambdas_vecreal = [];
plotting_bool = false;
mistakeThresh_real = 1E-3;
numSamples_int = 2000;
perturbAmount_vec = [ 0.01:0.02:0.2 0.3 0.5 ];
numsamples_vec = [1 5 10 50 100 500 1000 5000 ];
mistakeThresh_real_vec = [ 1E-3 5E-3 2E-2 2E-2 6E-2 1E-1 2E-1 4E-1 9E-1 2];
numMistakes_inv = nan( length(perturbAmount_vec), length(numsamples_vec), length(mistakeThresh_real_vec) );

for i=1:length(perturbAmount_vec)
    perturbAmount_real =  perturbAmount_vec(i);
    for j=1:length(numsamples_vec)
        numSamples_int = numsamples_vec(j);
        [ vecMistakes_avg_inv, vecMistakes_avg_lasso ] = script_univ( n_int, linear_bool, useLassoNmat_bool, Alimiter_real, numMatSamples_int,...
            perturbAmount_real, numSamples_int, mistakeThresh_real_vec, sigmas_vecreal,...
            lambdas_vecreal, plotting_bool );
        numMistakes_inv(i,j,:) = vecMistakes_avg_inv;
    end
end

% dimension1=perturbamount
% dimension2=numsamples_amount
% dimension3=mistakethresh_count

save('ws_data_updated')

