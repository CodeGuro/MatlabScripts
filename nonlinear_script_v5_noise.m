% %useful function: combnk(x,n) makes all combinations of n components from a
% vector 'x'

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
        product_den = 1; % product in denominator
        % size test (if the power set size is zero
        if size( vecJ, 2 ) < 1 
           % there are no nonzero elements for the current row
           new_vecX( i ) = alpha_null( i );
           vecX( i ) = alpha_null( i );
           continue; 
        end

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

    % matVecXItr (change in X with every iteration)
    for vecXItr = 1 : size( vecX, 1 )
        matVecXItr( vecXItr, iteration + 1 ) = vecX( vecXItr );
    end

    if iteration >= maxIterations
        break;
    end
    iteration = iteration + 1;
end

% we've found a solution
steady_vecX = vecX;
for steadyItr = 1 : size( steady_vecX, 1 )
    matSteadyVecX( :, recompute_itr + 1 ) = steady_vecX;
end

% display

for vecXItr = 1 : size( matVecXItr, 1 )
    plot( 1:size(matVecXItr,2), matVecXItr( vecXItr, : ) );
    xlabel( strcat( 'iteration, sigma=', num2str( sigma ) ) );
    ylabel( strcat( 'x_{', int2str( vecXItr ), '}' ) );
    pause(0.5);
end
clear matVecXItr;
    

for vecXItr = 1 : size( matSteadyVecX, 1 )
    plot( 1:size( matSteadyVecX, 2 ), matSteadyVecX( vecXItr, : ) );
    xlabel( 'sigma factor' );
    ylabel( strcat( 'x_', int2str( vecXItr ) ) );
    pause( 2 );
end