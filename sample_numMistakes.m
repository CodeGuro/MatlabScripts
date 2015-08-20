function [ vecMistakes_inv_avg, vecMistakes_inv_dev, vecMistakes_lasso_avg, vecMistakes_lasso_dev ] = sample_numMistakes( n, steady_vecX, matA, vecB, matK, matN, setsJ, powerSetsJ, ...
                alphas, alpha_null, itDiff_threshold, mistake_threshold, maxIterations, perturb_amount, num_samples, ...
                lambdas, sigma, linear, use_lasso_Nmat )

    vecMistakes_inv = nan( 1, num_samples );
    vecMistakes_lasso = nan( num_samples, length(lambdas) );

    for sample_num = 1:num_samples

        disp(['iterating... sample ' num2str(sample_num) ' out of ' num2str(num_samples)]);

        matdeltaX = construct_matdeltaX( n, steady_vecX, vecB, matK, matN, setsJ, powerSetsJ, ...
            alphas, alpha_null, itDiff_threshold, maxIterations, perturb_amount, sigma, linear );

        % We can now attempt to construct the linear matrix using the offsets
        matK_rec = matK_rec_useInv( n, matdeltaX );
        vecMistakes_inv( sample_num ) = nnz( logical( matK ) - logical( abs(matK_rec) > mistake_threshold ) );

        if length(lambdas) > 0
            matK_recs_lasso = matK_rec_useLasso( n, matdeltaX, lambdas, use_lasso_Nmat );
            for z=1:length(lambdas)
                numMistakes_lasso = nnz( logical( matA ) - logical( abs( matK_recs_lasso(:,:,z) ) > mistake_threshold ) );
                vecMistakes_lasso( sample_num, z ) = numMistakes_lasso;
            end
        end

    end
    
    % vecMistakes_inv is size: 1 x num_samples, vector of samples for a sigma
    vecMistakes_inv_avg = mean( vecMistakes_inv );
    vecMistakes_inv_dev = std( vecMistakes_inv );
    
    % vecMistakes_lasso is size: num_samples x length(lambdas), for a given sigma
    vecMistakes_lasso_avg = mean( vecMistakes_lasso, 1 );
    vecMistakes_lasso_dev = std( vecMistakes_lasso, 1 );
    

end
