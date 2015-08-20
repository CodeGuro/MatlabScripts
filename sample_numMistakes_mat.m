function [ vecMistakes_avg_allsigmas, vecMistakes_dev_allsigmas ] = sample_numMistakes_mat( n,...
    steady_vecX, matA, vecB, matK, matN, setsJ, powerSetsJ, alphas, ...
    alpha_null, itDiff_threshold, mistake_threshold, maxIterations, ...
    perturb_amount, num_samples, lambdas, sigmas, linear, use_lasso_Nmat )

    vecMistakes_avg_allsigmas = nan( length(sigmas), 1 + length(lambdas) );
    vecMistakes_dev_allsigmas = nan( length(sigmas), 1 + length(lambdas) );

    for sigma_it = 1:length(sigmas)
        disp(['sampling sigma: ' num2str(sigma_it) ' out of ' num2str(length(sigmas))]);
        
        sigma = sigmas( sigma_it );
        
        [ vecMistakes_avg, vecMistakes_dev ] ...
            = sample_numMistakes_sigma( n, steady_vecX, matA, vecB, matK, matN, setsJ, powerSetsJ, alphas, ...
            alpha_null, itDiff_threshold, mistake_threshold, maxIterations, perturb_amount, num_samples, ...
            lambdas, sigma, linear, use_lasso_Nmat );
        
        vecMistakes_avg_allsigmas( sigma_it, : ) = vecMistakes_avg;
        vecMistakes_dev_allsigmas( sigma_it, : ) = vecMistakes_dev;
        
    end

end

