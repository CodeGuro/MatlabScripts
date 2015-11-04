linear_bool = false;
useLasso_bool = true;
useQLasso_bool = true;
n_int = 10;
Alimiter_real = 0.9;
numMatSamplesInternal_int = 1;
numMatSamplesExternal_int = 10;
sigmas_vecreal = [ 1E-2 ];
lambdas_vecreal = [ 0 1E-5:1E-5:1E-4 2E-4:1E-4:1E-3 2E-3:1E-3:1E-2 2E-2:1E-2:1E-1 ];
plotting_bool = false;
perturbAmount_vec = [ 1E-3:1E-3:1E-2 2E-2:1E-2:1E-1 2E-1:1E-1:1E-0 ];
numsamples_vec = [ 1:9 10:10:100 200:100:1000 2000:1000:10000 20000:10000:100000 ];
mistakeThresh_real_vec = [ 1E-6:1E-6:1E-5 2E-5:1E-5:1E-4 2E-4:1E-4:1E-3 2E-3:1E-3:1E-2 2E-2:1E-2:1E-1 2E-1:1E-1:1E-0];
numMistakes_inv = nan( length(perturbAmount_vec), length(numsamples_vec), length(mistakeThresh_real_vec), numMatSamplesExternal_int );
numMistakes_lasso = nan( length(perturbAmount_vec), length(numsamples_vec), length(lambdas_vecreal), numMatSamplesExternal_int );
numMistakes_Qlasso = numMistakes_lasso;
nonoise_inv = nan( length( perturbAmount_vec ), length( mistakeThresh_real_vec ), numMatSamplesExternal_int );
nonoise_lasso = nan( length( perturbAmount_vec ), length( lambdas_vecreal ), numMatSamplesExternal_int ); 
nonoise_Qlasso = nonoise_lasso;

for matsample = 1:numMatSamplesExternal_int
    disp( ['sampling matrix ' num2str(matsample) ' of ' num2str(numMatSamplesExternal_int) ] );
    Structargs = gen_vars( n_int, Alimiter_real );
    for i=1:length(perturbAmount_vec)
        perturbAmount_real =  perturbAmount_vec(i);
        [nonoise_inv_vec, nonoise_lasso_vec, nonoise_Qlasso_vec ] = script_univ( n_int, linear_bool, useLasso_bool, useQLasso_bool, Alimiter_real, numMatSamplesInternal_int,...
            perturbAmount_real, 1, mistakeThresh_real_vec, 0,...
            lambdas_vecreal, plotting_bool, Structargs );
        nonoise_inv( i, :, matsample ) = nonoise_inv_vec;
        nonoise_lasso( i, :, matsample ) = nonoise_lasso_vec;
        nonoise_Qlasso( i, :, matsample ) = nonoise_Qlasso_vec;
        disp( [ 'sampling perturbation ' num2str(i) ' of ' num2str(length(perturbAmount_vec)) ': ' num2str( perturbAmount_real ) ] );
        for j=1:length(numsamples_vec)
            numSamples_int = numsamples_vec(j);
            disp( ['sampling repitition ' num2str(j) ' of ' num2str(length(numsamples_vec)) ': ' num2str( numSamples_int ) ]);
            [ vecMistakes_avg_inv, vecMistakes_avg_lasso, vecMistakes_avg_Qlasso ] = script_univ( n_int, linear_bool, useLasso_bool, useQLasso_bool, Alimiter_real, numMatSamplesInternal_int,...
                perturbAmount_real, numSamples_int, mistakeThresh_real_vec, sigmas_vecreal,...
                lambdas_vecreal, plotting_bool, Structargs );
            numMistakes_inv(i,j,:,matsample) = vecMistakes_avg_inv;
            numMistakes_lasso(i,j,:,matsample) = vecMistakes_avg_lasso;
            numMistakes_Qlasso(i,j,:,matsample) = vecMistakes_avg_Qlasso;
        end
    end
end

numMistakes_inv_avg = mean( numMistakes_inv, 4 );
numMistakes_lasso_avg = mean( numMistakes_lasso, 4 );
numMistakes_Qlasso_avg = mean( numMistakes_Qlasso, 4 );

save('ws_data_3.mat');
% dimension1=perturbamount
% dimension2=numsamples_amount
% dimension3=mistakethresh_count
% fname_lasso = 'lasso';
% fname_inv = 'linear inverse';
% mkdir( fname_lasso );
% mkdir( fname_inv );

