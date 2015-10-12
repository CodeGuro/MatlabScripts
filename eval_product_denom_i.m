function product_denom_i = eval_product_denom_i( i, vecX, matK, matN, setsJ )
    product_denom_i = 1; % default
        for j = 1 : size( setsJ{ i }, 2 )
            jval = setsJ{ i }{ j };
            product_denom_i = product_denom_i * ( 1 + epsilonfunc( i, jval, vecX, matK, matN ) );
        end
end