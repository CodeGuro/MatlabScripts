function [ setsJ, powerSetsJ, alphas, alpha_null, matK, matN ] = gen_nonlin_vars( n, matA )

    matN = matA * 2;
    matK = zeros( n, n );
    alpha_null = rand( n, 1 );
    
    for i = 1 : n
        vecIdx = find( matA( i, : ) );
        pSet = {};
        setsJ{ i } = {};
        for j = 1 : length( vecIdx )
            % fill the K(i,j) uniformly [0,1] where matA(i,j)==1
            matK( i, vecIdx( j ) ) = rand();
            % for each row in matA, define J = { j | matA(i,j)==1 }
            % and P(J)
            setsJ{ i }{ j } = vecIdx( j );
            s = vecsToSets( combnk( vecIdx, j ) );
            pSet = appendToSet( pSet, s );
        end
        powerSetsJ{ i } = pSet;
    end

    % generate alphas
    alphas = cell( size( powerSetsJ ) );
    for i = 1 : size( powerSetsJ, 2 )
        for j = 1 : size( powerSetsJ{ i }, 2 )
            alphas{ i }{ j } = rand();
        end
    end

end

