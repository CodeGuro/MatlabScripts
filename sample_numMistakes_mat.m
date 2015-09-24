function [ vecMistakes_avg_nsamples, vecMistakes_dev_nsamples ] = sample_numMistakes_mat( n,...
    steady_vecX, matA, vecB, matK, matN, setsJ, powerSetsJ, alphas, ...
    alpha_null, itDiff_threshold, mistake_threshold, maxIterations, ...
    perturb_amount, num_samples, lambdas, sigmas, linear, use_lasso_Nmat )

    vecMistakes_avg_nsamples = nan( length(sigmas), length(mistake_threshold) );
    vecMistakes_dev_nsamples = nan( length(sigmas), length(mistake_threshold) );

    for sigma_it = 1:length(sigmas)
        disp(['sampling sigma: ' num2str(sigma_it) ' out of ' num2str(length(sigmas))]);
        
        sigma = sigmas( sigma_it );
        
        [ vecMistakes_avg, vecMistakes_dev ] ...
            = sample_numMistakes_sigma( n, steady_vecX, matA, vecB, matK, matN, setsJ, powerSetsJ, alphas, ...
            alpha_null, itDiff_threshold, mistake_threshold, maxIterations, perturb_amount, num_samples, ...
            lambdas, sigma, linear, use_lasso_Nmat );
        
        vecMistakes_avg_nsamples( sigma_it, : ) = vecMistakes_avg;
        vecMistakes_dev_nsamples( sigma_it, : ) = vecMistakes_dev;
        
    end

end

