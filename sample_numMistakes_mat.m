function [ vecMistakes_inv_nsamples, vecMistakes_lasso_nsamples, vecMistakes_Qlasso_nsamples ] = sample_numMistakes_mat( n,...
    steady_vecX, matA, vecB, matK, matN, setsJ, powerSetsJ, alphas, ...
    alpha_null, itDiff_threshold, mistake_threshold, maxIterations, ...
    perturb_amount, num_samples, lambdas, sigmas, linear, use_lasso, use_Qlasso )

    vecMistakes_inv_nsamples = nan( length(sigmas), length(mistake_threshold) );
    vecMistakes_lasso_nsamples = nan( length(sigmas), length(lambdas) );
    vecMistakes_Qlasso_nsamples = vecMistakes_lasso_nsamples;

    for sigma_it = 1:length(sigmas)
        disp(['sampling sigma: ' num2str(sigma_it) ' out of ' num2str(length(sigmas))]);
        
        sigma = sigmas( sigma_it );
        
        [ vecMistakes_inv, vecMistakes_lasso, vecMistakes_Qlasso ] ...
            = sample_numMistakes_sigma( n, steady_vecX, matA, vecB, matK, matN, setsJ, powerSetsJ, alphas, ...
            alpha_null, itDiff_threshold, mistake_threshold, maxIterations, perturb_amount, num_samples, ...
            lambdas, sigma, linear, use_lasso, use_Qlasso );
        
        vecMistakes_inv_nsamples( sigma_it, : ) = vecMistakes_inv;
        vecMistakes_lasso_nsamples( sigma_it, : ) = vecMistakes_lasso;
        vecMistakes_Qlasso_nsamples( sigma_it, : ) = vecMistakes_Qlasso;
        
    end

end

