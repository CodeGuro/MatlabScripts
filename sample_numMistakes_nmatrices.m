function [ vecMistakes_avg_inv, vecMistakes_avg_lasso ] = sample_numMistakes_nmatrices( n, ...
            itDiff_threshold, mistake_threshold, maxIterations, numMatrix_samples, ...
            perturb_amount, num_samples, lambdas, sigmas, linear, use_lasso_Nmat, ...
            A_limiter, Structargs )

    matMistakes_inv_sigmaxthresh = nan( length(sigmas), length(mistake_threshold), numMatrix_samples );
    matMistakes_lasso_sigmaxlambdas = nan( length(sigmas), length(lambdas), numMatrix_samples );
    for M_sample = 1:numMatrix_samples

        disp( [ 'sampling matrix: ' num2str(M_sample) ' out of ' num2str(numMatrix_samples) ] );

        % generate setsJ, powerSetsJ, and matK, and other nonlinear vars
        [ matA, vecB, setsJ, powerSetsJ, alphas, alpha_null, matK, matN ] = gen_vars( n, A_limiter );

        % find the steady state
        if linear
            steady_vecX = findSteady_lin( n, matK, vecB );
        else
            steady_vecX = findSteady_nonlin( n, matK, matN, setsJ, powerSetsJ, alphas, ...
                alpha_null, itDiff_threshold, maxIterations );
        end

        [ vecMistakes_inv_nsamples, vecMistakes_lasso_nsamples ] = sample_numMistakes_mat( n,...
            steady_vecX, matA, vecB, matK, matN, setsJ, powerSetsJ, alphas, ...
            alpha_null, itDiff_threshold, mistake_threshold, maxIterations, ...
            perturb_amount, num_samples, lambdas, sigmas, linear, use_lasso_Nmat );

        matMistakes_inv_sigmaxthresh( :, :, M_sample ) = vecMistakes_inv_nsamples;
        matMistakes_lasso_sigmaxlambdas( :, :, M_sample ) = vecMistakes_lasso_nsamples;

    end
    
    % dimension 1: sigmas
    % dimension 2: mistake_threshold || lambdas 
    % dimension 3: matrix samples
    vecMistakes_avg_inv = mean( matMistakes_inv_sigmaxthresh, 3 );
    vecMistakes_avg_lasso = mean( matMistakes_lasso_sigmaxlambdas, 3 );

end

