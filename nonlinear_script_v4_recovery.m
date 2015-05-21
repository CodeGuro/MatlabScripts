% generate the random n x n matrix
n = 3;
matA = [0 1 0; 0 0 1; 1 0 0]; %randi( [0;1], n, n ); %[0,1,1;1,0,0;1,1,0];%
matA( logical( eye( n ) ) ) = 0;
matK = zeros( n, n );
matN = matA * 2;
vecX = zeros( n, 1 ); % the gene expressions
new_vecX = zeros( n, 1 ); % refresh
setsJ = {};
powerSetsJ = {};
alphas = {};

for i = 1 : n
    vecIdx = find( matA( i, : ) );
    pSet = {};
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
for i = 1 : size( powerSetsJ, 2 )
    for j = 1 : size( powerSetsJ{ i }, 2 )
        alphas{ i }{ j } = rand();
    end
end

alpha_null = rand( n, 1 );

    
vecX = zeros( n, 1 );
new_vecX = zeros( n, 1 );

iteration = 1;
maxIterations = 1000;
finished = false;
while ~finished

    for i = 1 : n
        vecJ = find( matA( i, : ) );
        sum_num = 0; % sum in the numerator

        % begin evaluating the sum (numerator)
        for j = 1 : size( powerSetsJ{ i }, 2 )
            product = 1;
            for it = 1 : size( powerSetsJ{ i }{ j }, 2 )
                vec = powerSetsJ{ i }{ j };
                alpha = alphas{ i }{ j };
                for it2 = 1 : size( vec, 2 )
                    product = product * epsilonfunc( i, vec( it2 ), vecX, matK, matN );
                end
                product = product * alpha;
            end
            sum_num = sum_num + product;
        end

        sum_num = sum_num + alpha_null( i );

        %now to begin the denominator product
        product_den = 1; % default
        for j = 1 : size( setsJ{ i }, 2 )
            jval = setsJ{ i }{ j };
            product_den = product_den * ( 1 + epsilonfunc( i, vec( it2 ), vecX, matK, matN ) );
        end

        new_vecX( i ) = sum_num / product_den;

    end

    if size( find( abs( vecX - new_vecX ) > 0.001 ), 1 ) == 0
        finished=true;
        vecX-new_vecX %display difference
    end

    vecX = new_vecX;

    if iteration >= maxIterations
        break;
    end
    iteration = iteration + 1;
end

% now that we've found a solution, we can start the perturbations

perturb_amount = -0.1;
steady_vecX = vecX;
matdeltaX = NaN( n, n );


for p = 1:n % perturbation index
    finished = false;
    
    new_vecX = zeros( n, 1 );
    vecX= zeros( n, 1 );
    new_vecX( p ) = steady_vecX( p ) + perturb_amount;
    vecX( p ) = steady_vecX( p ) + perturb_amount;
    
    while ~finished

        for i = 1 : n
            if i ~= p
                
                vecJ = find( matA( i, : ) );
                sum_num = 0; % sum in the numerator

                % begin evaluating the sum (numerator)
                for j = 1 : size( powerSetsJ{ i }, 2 )
                    product = 1;
                    for it = 1 : size( powerSetsJ{ i }{ j }, 2 )
                        vec = powerSetsJ{ i }{ j };
                        alpha = alphas{ i }{ j };
                        % iterator through the vector (i.e. [1,3] from example)
                        for it2 = 1 : size( vec, 2 )
                            product = product * epsilonfunc( i, vec( it2 ), vecX, matK, matN );
                        end
                        product = product * alpha;
                    end
                    sum_num = sum_num + product;
                end

                sum_num = sum_num + alpha_null( i );
                %sum_num is now the complete for row i

                %now to begin the denominator product
                product_den = 1; % default
                for j = 1 : size( setsJ{ i }, 2 )
                    jval = setsJ{ i }{ j };
                    product_den = product_den * ( 1 + epsilonfunc( i, vec( it2 ), vecX, matK, matN ) );
                end

                new_vecX( i ) = sum_num / product_den;
            end
        end

        if size( find( abs( vecX - new_vecX ) > 1E-3 ), 1 ) == 0
            finished=true;
            vecDiff = vecX-new_vecX %steady state should differ by below 1E-3
            
            vecdeltaX = new_vecX - steady_vecX;
            matdeltaX( p, : ) = vecdeltaX;
        end

    vecX = new_vecX;

    end
    
end


% We can now attempt to construct the linear matrix using the offsets
% (deltamatX)
lin_mat = nan( n, n );
for current = 1 : n
    selection = setdiff( 1:n, current );
    lin_mat_cur_row = matdeltaX( selection, selection ) \ matdeltaX( selection, current );
    lin_mat_cur_row = insert( lin_mat_cur_row, 0, current );
    lin_mat( current, : ) = lin_mat_cur_row;
end

matA_recovered = abs(lin_mat) > 1E-2;
matA_logical = logical( matA );
