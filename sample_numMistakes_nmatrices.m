function [ numMistakes_avg_Msamples, numMistakes_devs_Msamples ] = sample_numMistakes_nmatrices( n, ...
            itDiff_threshold, mistake_threshold, maxIterations, numMatrix_samples, ...
            perturb_amount, num_samples, lambdas, sigmas, linear, use_lasso_Nmat, ...
            A_limiter )

    vecMistakes_avg_nsamples = nan( length(sigmas), 1 + length(lambdas), numMatrix_samples );
    vecMistakes_dev_nsamples = nan( length(sigmas), 1 + length(lambdas), numMatrix_samples );
    for M_sample = 1:numMatrix_samples

        disp( [ 'sampling matrix: ' num2str(M_sample) ' out of ' num2str(numMatrix_samples) ] );

        matA = rand(n,n) > A_limiter;
        matA( logical( eye( n ) ) ) = 0;
        b_min = 1E-5;
        b_max = 1;
        vecB = b_min + (b_max - b_min).*rand( [n,1] );
        % generate setsJ, powerSetsJ, and matK, and other nonlinear vars
        [ setsJ, powerSetsJ, alphas, alpha_null, matK, matN ] = gen_nonlin_vars( n, matA );

        % find the steady state
        if linear
            steady_vecX = findSteady_lin( n, matK, vecB );
        else
            steady_vecX = findSteady_nonlin( n, matK, matN, setsJ, powerSetsJ, alphas, ...
                alpha_null, itDiff_threshold, maxIterations );
        end

        [ vecMistakes_avg_allsigmas, vecMistakes_dev_allsigmas ] = sample_numMistakes_mat( n,...
            steady_vecX, matA, vecB, matK, matN, setsJ, powerSetsJ, alphas, ...
            alpha_null, itDiff_threshold, mistake_threshold, maxIterations, ...
            perturb_amount, num_samples, lambdas, sigmas, linear, use_lasso_Nmat );

        vecMistakes_avg_nsamples( :, :, M_sample ) = vecMistakes_avg_allsigmas;
        vecMistakes_dev_nsamples( :, :, M_sample ) = vecMistakes_dev_allsigmas;

    end
    
    % dimension 1: sigmas
    % dimension 2: [ inv, lambdas ]
    numMistakes_avg_Msamples = mean( vecMistakes_avg_nsamples, 3 );
    numMistakes_devs_Msamples = mean( vecMistakes_dev_nsamples, 3 );

end

