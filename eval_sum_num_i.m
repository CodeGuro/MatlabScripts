function sum_num_i = eval_sum_num_i( i, vecX, matK, matN, powerSetsJ, alphas, alpha_null )
    sum_num_i = 0; % sum in the numerator
    % begin evaluating the sum (numerator)
    for j = 1 : size( powerSetsJ{ i }, 2 )
        product = 1;
        vec = powerSetsJ{ i }{ j };
        alpha = alphas{ i }{ j };
        for it2 = 1 : size( vec, 2 )
            product = product * epsilonfunc( i, vec( it2 ), vecX, matK, matN );
        end
        product = product * alpha;
        sum_num_i = sum_num_i + product;
    end
    sum_num_i = sum_num_i + alpha_null( i );
end

