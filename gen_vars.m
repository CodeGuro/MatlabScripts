function [ Structargs ] = gen_vars( n, A_limiter )

    Structargs.matA = rand(n,n) > A_limiter;
    Structargs.matA( logical( eye( n ) ) ) = 0;
    b_min = 1E-5;
    b_max = 1;
    Structargs.vecB = b_min + (b_max - b_min).*rand( [n,1] );

    Structargs.matN = Structargs.matA * 2;
    Structargs.matK = zeros( n, n );
    Structargs.alpha_null = rand( n, 1 );
    
    for i = 1 : n
        vecIdx = find( Structargs.matA( i, : ) );
        pSet = {};
        Structargs.setsJ{ i } = {};
        for j = 1 : length( vecIdx )
            % fill the K(i,j) uniformly [0,1] where matA(i,j)==1
            Structargs.matK( i, vecIdx( j ) ) = rand();
            % for each row in matA, define J = { j | matA(i,j)==1 }
            % and P(J)
            Structargs.setsJ{ i }{ j } = vecIdx( j );
            s = vecsToSets( combnk( vecIdx, j ) );
            pSet = appendToSet( pSet, s );
        end
        Structargs.powerSetsJ{ i } = pSet;
    end

    % generate alphas
    Structargs.alphas = cell( size( Structargs.powerSetsJ ) );
    for i = 1 : size( Structargs.powerSetsJ, 2 )
        for j = 1 : size( Structargs.powerSetsJ{ i }, 2 )
            Structargs.alphas{ i }{ j } = rand();
        end
    end    
    
end

