function [ matdeltaX ] = construct_matdeltaX( n, steady_vecX, vecB, matK, matN, setsJ, powerSetsJ, ...
    alphas, alpha_null, itDiff_threshold, maxIterations, perturb_amount, sigmas, sigma_it, linear )

    matdeltaX = NaN( n, n );
    for p = 1:n % perturbation index

        disp(['perturbation index: ' num2str(p) ' out of ' num2str(n)]);

        if linear
            vecX_P = perturb_lin_v2( n, steady_vecX, matK, vecB, p, perturb_amount, sigmas( sigma_it ) );
        else
            vecX_P = perturb_nonlin( n, steady_vecX, matK, matN, setsJ, powerSetsJ, alphas, ...
                alpha_null, itDiff_threshold, maxIterations, p, perturb_amount, sigmas( sigma_it ) );
        end

        matdeltaX( :, p ) = vecX_P - steady_vecX;
    end
    
end

