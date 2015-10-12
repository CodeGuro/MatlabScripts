function vecX_P = perturb_nonlin( n, steady_vecX, matK, matN, setsJ, powerSetsJ, alphas, alpha_null, itDiff_threshold, maxIterations, p, perturb_amount, sigma )

    vecX = zeros( n, 1 );
    finished = false;
    iteration = 1;
    while ~finished

        next_vecX = nonlinear_func( n, vecX, matK, matN, setsJ, powerSetsJ, alphas, alpha_null );
        next_vecX( p ) = steady_vecX( p ) + perturb_amount;

        if size( find( abs( vecX - next_vecX ) > itDiff_threshold ), 1 ) == 0 || iteration >= maxIterations
            finished=true;
            if iteration >= maxIterations
                disp('warning! max iteration exceeded! aborting...');
            end
            vecDiff = vecX-next_vecX;

            noise = randn(n, 1) * sigma;
            noise(p) = 0;
            perturbed_i_steady_vecX = next_vecX;
            perturbed_i_steady_vecX( p ) = steady_vecX( p ) + perturb_amount;
            perturbed_i_steady_vecX = perturbed_i_steady_vecX + noise; % noise added here
        else
            vecX = next_vecX;
        end
        iteration = iteration + 1;
    end
    
    vecX_P = perturbed_i_steady_vecX;
    
end

