function steady_vecX = findSteady_nonlin( n, matK, matN, setsJ, powerSetsJ, alphas, alpha_null, itDiff_threshold, maxIterations )

    vecX = zeros( n, 1 );
    finished = false;
    iteration = 1;
    while ~finished
        next_vecX = nonlinear_func( n, vecX, matK, matN, setsJ, powerSetsJ, alphas, alpha_null );

        if size( find( abs( vecX - next_vecX ) > itDiff_threshold ), 1 ) == 0 || iteration >= maxIterations
            finished=true;
            if iteration >= maxIterations
                disp('warning! max iteration exceeded!');
            end
            vecDiff = vecX-next_vecX;
        else
            vecX = next_vecX;
        end

        iteration = iteration + 1;
    end
    
    steady_vecX = next_vecX;
end

