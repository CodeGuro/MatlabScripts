% options
linear = false;
use_lasso_Nmat = false;
n = 5;
A_limiter = 0.8;
numMatrix_samples = 6;
maxIterations = 1000;
itDiff_threshold = 1E-4;
perturb_amount = 1;
num_samples = 20;
mistake_threshold = 1E-1;
sigmas = 0:0.03:0.5;
lambdas = [2E-4 2E-2 0.1 0.2 1 2E10 ];


% algorithm starts here
for M_sample = 1:numMatrix_samples

    matA = rand(n,n) > A_limiter;
    matA( logical( eye( n ) ) ) = 0;
    b_min = 1E-5;
    b_max = 1;
    vecB = b_min + (b_max - b_min).*rand( [n,1] );
    
    if linear
        func_type = 'linear';
    else
        func_type = 'nonlinear';
    end

    % generate setsJ, powerSetsJ, and matK, and other nonlinear vars
    [ setsJ, powerSetsJ, alphas, alpha_null, matK, matN ] = gen_nonlin_vars( n, matA );

    % find the steady state
    if linear
        steady_vecX = findSteady_lin( n, matK, vecB );
    else
        steady_vecX = findSteady_nonlin( n, matK, matN, setsJ, powerSetsJ, alphas, ...
            alpha_null, itDiff_threshold, maxIterations );
    end

    % now that we've found a steady state, we can start the perturbations
    
    for sigma_it = 1:length(sigmas)
        disp(['sampling sigma: ' num2str(sigma_it) ' out of ' num2str(length(sigmas))]);
        
        sigma = sigmas( sigma_it );
        
        [ vecMistakes_inv_avg, vecMistakes_inv_dev, vecMistakes_lasso_avg, vecMistakes_lasso_dev ] ...
            = sample_numMistakes( n, steady_vecX, matA, vecB, matK, matN, setsJ, powerSetsJ, alphas, ...
            alpha_null, itDiff_threshold, mistake_threshold, maxIterations, perturb_amount, num_samples, ...
            lambdas, sigma, linear, use_lasso_Nmat );
        
        vecMistakes_inv_msamples( M_sample, sigma_it ) = vecMistakes_inv_avg; % std(dim1) on assgined var for samples across matrices
        vecstDevs_inv_msamples( M_sample, sigma_it ) = vecMistakes_inv_dev; % std(dim1) on assigned var for samples across matrices

        if length(lambdas) > 0
            vecMistakes_lasso_msamples( :, sigma_it, M_sample ) = vecMistakes_lasso_avg; % std(dim3) on assigned var for samples across matrices
            vecstDevs_lasso_msamples( :, sigma_it, M_sample ) = vecMistakes_lasso_dev; % std(dim3) on assigned var for samples across matrices
        end

    end
    
end

vecMistakes_inv_msamples_avg = mean( vecMistakes_inv_msamples, 1 );
vecstDevs_inv_msamples_avg = mean( vecstDevs_inv_msamples, 1 );
vecMistakes_avg_lasso = mean( vecMistakes_lasso_msamples, 3 );
vecstDevs_lasso = mean( vecstDevs_lasso_msamples, 3 );

% output
plot_x = (ones(length(lambdas),1)*sigmas)';
plot_y = vecMistakes_avg_lasso';
plot_devs = vecstDevs_lasso';
legends = cell( 1, length(lambdas) );
for lambda_it = 1:length(lambdas)
    legends{ lambda_it } = ['lasso: lambda=' num2str( lambdas(lambda_it) ) ];
end

plot_x( :, end+1 ) = sigmas;
plot_y(:, end+1 ) = vecMistakes_inv_msamples_avg;
plot_devs(:, end+1 ) = vecstDevs_inv_msamples_avg;
legends{length(lambdas) + 1} = 'inv';

errorbar( plot_x, plot_y, plot_devs, ':o' );
legend( 'avg mistakes' );
xlabel('sigma');
ylabel('recovery mistakes');
legend( legends );
title( {[num2str(n) 'X' num2str(n) ' matrix of a ' func_type ' regulatory network, with '...
    num2str(perturb_samples) ' samples per perturbation'],...
    ['perturb amount=' num2str(perturb_amount) ', '...
    'A limiter=' num2str(A_limiter) ', ' 'nnz(matA)=' num2str(nnz(matA)) ', and ' num2str(numMatrix_samples) ' matrices sampled']} );
%axis( [ sigmas(1) sigmas(end) 0 10] );
axis('auto');
