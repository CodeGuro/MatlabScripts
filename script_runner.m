linear_bool = false;
useLassoNmat_bool =  false;
n_int = 10;
Alimiter_real = 0.9;
numMatSamplesInternal_int = 1;
numMatSamplesExternal_int = 30;
sigmas_vecreal = [ 0.04 ];
lambdas_vecreal = 0:0.001:0.1;
plotting_bool = false;
mistakeThresh_real = 1E-3;
numSamples_int = 2000;
perturbAmount_vec = [ 0 3E-4:2E-3:0.01 0.01:0.02:0.2 0.3 0.5 ];
numsamples_vec = [1 5 10 50 100 500 1000 3000 5000 ];
mistakeThresh_real_vec = [ 0:2E-3:0.2 2];
numMistakes_inv = nan( length(perturbAmount_vec), length(numsamples_vec), length(mistakeThresh_real_vec), numMatSamplesExternal_int );
numMistakes_lasso = nan( length(perturbAmount_vec), length(numsamples_vec), length(lambdas_vecreal), numMatSamplesExternal_int );

for matsample = 1:numMatSamplesExternal_int
    disp( ['sampling matrix ' num2str(matsample) ' of ' num2str(numMatSamplesExternal_int) ] );
    Structargs = gen_vars( n_int, Alimiter_real );
    for i=1:length(perturbAmount_vec)
        perturbAmount_real =  perturbAmount_vec(i);
        disp( [ 'sampling perturbation ' num2str(i) ' of ' num2str(length(perturbAmount_vec)) ': ' num2str( perturbAmount_real ) ] );
        for j=1:length(numsamples_vec)
            numSamples_int = numsamples_vec(j);
            disp( ['sampling repitition ' num2str(j) ' of ' num2str(length(numsamples_vec)) ': ' num2str( numSamples_int ) ]);
            [ vecMistakes_avg_inv, vecMistakes_avg_lasso ] = script_univ( n_int, linear_bool, useLassoNmat_bool, Alimiter_real, numMatSamplesInternal_int,...
                perturbAmount_real, numSamples_int, mistakeThresh_real_vec, sigmas_vecreal,...
                lambdas_vecreal, plotting_bool, Structargs );
            numMistakes_inv(i,j,:,matsample) = vecMistakes_avg_inv;
            numMistakes_lasso(i,j,:,matsample) = vecMistakes_avg_lasso;
        end
    end
end

numMistakes_inv_avg = mean( numMistakes_inv, 4 );
numMistakes_lasso_avg = mean( numMistakes_lasso, 4 );

% dimension1=perturbamount
% dimension2=numsamples_amount
% dimension3=mistakethresh_count

save('ws_data_updated')

