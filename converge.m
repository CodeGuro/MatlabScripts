function next_vecX = converge( n, vecX, matK, matN, setsJ, powerSetsJ, alphas, alpha_null )
    next_vecX = zeros( n, 1 );
    for i = 1 : n
        sum_num = eval_sum_num_i( i, vecX, matK, matN, powerSetsJ, alphas, alpha_null );
        product_den = eval_product_denom_i( i, vecX, matK, matN, setsJ );
        next_vecX( i ) = sum_num / product_den;
    end
end

