function [ vecMistakes_avg_inv, vecMistakes_avg_lasso, vecMistakes_avg_Qlasso ] = sample_numMistakes_nmatrices( n, ...
            itDiff_threshold, mistake_threshold, maxIterations, numMatrix_samples, ...
            perturb_amount, num_samples, lambdas, sigmas, linear, use_lasso, use_Qlasso, ...
            A_limiter, Structargs_param )

    matMistakes_inv_sigmaxthresh = nan( length(sigmas), length(mistake_threshold), numMatrix_samples );
    matMistakes_lasso_sigmaxlambdas = nan( length(sigmas), length(lambdas), numMatrix_samples );
    matMistakes_Qlasso_sigmaxlambdas = matMistakes_lasso_sigmaxlambdas;
    for M_sample = 1:numMatrix_samples

        % disp( [ 'sampling matrix: ' num2str(M_sample) ' out of ' num2str(numMatrix_samples) ] );

        % generate setsJ, powerSetsJ, and matK, and other nonlinear vars
        if isempty(fieldnames(Structargs_param))
            Structargs = gen_vars( n, A_limiter );
        else
            Structargs = Structargs_param;
        end
        fields = fieldnames( Structargs );
        for field_it = 1:length(fields)
            field = fields{field_it};
            eval(sprintf([field '=%s;'], 'Structargs.(field)'));
        end

        % find the steady state
        if linear
            steady_vecX = findSteady_lin( n, matK, vecB );
        else
            steady_vecX = findSteady_nonlin( n, matK, matN, setsJ, powerSetsJ, alphas, ...
                alpha_null, itDiff_threshold, maxIterations );
        end

        [ vecMistakes_inv_nsamples, vecMistakes_lasso_nsamples, ...
            vecMistakes_Qlasso_nsamples ] = sample_numMistakes_mat( n,...
            steady_vecX, matA, vecB, matK, matN, setsJ, powerSetsJ, alphas, ...
            alpha_null, itDiff_threshold, mistake_threshold, maxIterations, ...
            perturb_amount, num_samples, lambdas, sigmas, linear, use_lasso, use_Qlasso );

        matMistakes_inv_sigmaxthresh( :, :, M_sample ) = vecMistakes_inv_nsamples;
        matMistakes_lasso_sigmaxlambdas( :, :, M_sample ) = vecMistakes_lasso_nsamples;
        matMistakes_Qlasso_sigmaxlambdas( :, :, M_sample ) = vecMistakes_Qlasso_nsamples;

    end
    
    % dimension 1: sigmas
    % dimension 2: mistake_threshold || lambdas 
    % dimension 3: matrix samples
    vecMistakes_avg_inv = mean( matMistakes_inv_sigmaxthresh, 3 );
    vecMistakes_avg_lasso = mean( matMistakes_lasso_sigmaxlambdas, 3 );
    vecMistakes_avg_Qlasso = mean( matMistakes_Qlasso_sigmaxlambdas, 3 );

end

