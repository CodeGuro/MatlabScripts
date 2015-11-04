function [ vecMistakes_inv, vecMistakes_lasso, vecMistakes_Qlasso ] = sample_numMistakes_sigma( n, steady_vecX, matA, vecB, matK, matN, setsJ, powerSetsJ, ...
                alphas, alpha_null, itDiff_threshold, mistake_threshold, maxIterations, perturb_amount, num_samples, ...
                lambdas, sigma, linear, use_lasso, use_Qlasso )

    vecMistakes_inv = nan( 1, length(mistake_threshold) );
    vecMistakes_lasso = nan( 1, length(lambdas) );
    vecMistakes_Qlasso = vecMistakes_lasso;
    
    % noise added outside for speedy algorithm
    matdeltaX = construct_matdeltaX( n, steady_vecX, vecB, matK, matN, setsJ, powerSetsJ, ...
            alphas, alpha_null, itDiff_threshold, maxIterations, perturb_amount, 0, linear );
        
    noisy_matdeltaX_tiled = repmat( matdeltaX, 1,  num_samples ) + ( repmat(ones(n)-eye(n), 1, num_samples) .* (randn( n, n*num_samples) .* sigma) );
    
    matK_rec = matK_rec_useInv( n, noisy_matdeltaX_tiled, num_samples );
    
    for it=1:length(mistake_threshold)
        vecMistakes_inv( it ) = nnz( matA - logical( abs(matK_rec) > mistake_threshold(it) ) );
    end
    % disp( ['mistake count: ' num2str(vecMistakes_inv(1)) ] );
    
    % normal lasso
    if ~isempty(lambdas) && use_lasso
        matK_recs_lasso = matK_rec_useLasso( n, noisy_matdeltaX_tiled, num_samples, lambdas, false );
        for z=1:length(lambdas)
            numMistakes_lasso = nnz( logical( matA ) - logical( abs( matK_recs_lasso(:,:,z) ) ) );
            vecMistakes_lasso( 1, z ) = numMistakes_lasso;
        end
    end
    
    % quadratic lasso
    if ~isempty(lambdas) && use_Qlasso
        matK_recs_Qlasso = matK_rec_useLasso( n, noisy_matdeltaX_tiled, num_samples, lambdas, true );
        for z=1:length(lambdas)
            numMistakes_Qlasso = nnz( logical( matA ) - logical( abs( matK_recs_Qlasso(:,:,z) ) ) );
            vecMistakes_Qlasso( 1, z ) = numMistakes_Qlasso;
        end
    end

end

