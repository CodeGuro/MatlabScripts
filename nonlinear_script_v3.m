%%%%%   notes   %%%%%
% %useful function: combnk(x,n) makes all combinations of n components from a
% vector 'x'

% % the set
% x = 1:4;
% for n = 1:4
% combnz = combnk(x,n);
% %or% combnz{nn} = combnk(x,nn);
% end
%%%%%   %%%%%   %%%%%

% generate the random nxn matrix
n = 4;
matA = randi( [0;1], n, n ); %[0,1,1;1,0,0;1,1,0];%
matA( logical( eye( n ) ) ) = 0;
matK = zeros( n, n );
matN = matA * 2; % since n_ij=1 where matA_ij=1 (set to a lower value for higher size vector
vecX = zeros(n,1); % this will be refreshed with every iteration
new_vecX = zeros( n, 1 );
setsJ = {};
powerSetsJ = {};
alphas = {};
finished = false;

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
        % rown [ 0.4, 0, 0.8 ] ==> set {1,3,[1,3]}
        % alpha_rown,1*epsilon_rown,1 + alpha_rown,3 * alpha*epsilon_rown,3
        % + alpha_rown_[1,3] * alpha_epsilon_rown,1 * alpha_epsilon_rown,3
        for j = 1 : size( powerSetsJ{ i }, 2 )
            product = 1;
            for it = 1 : size( powerSetsJ{ i }{ j }, 2 )
                vec = powerSetsJ{ i }{ j };
                alpha = alphas{ i }{ j };
                % iterator through the vector (i.e. [1,3] from example)
                for it2 = 1 : size( vec, 2 )
                    product = product * power( vecX( vec( it2 ) ) / matK( vec( it2 ), i ), matN( vec( it2 ), i ) );
                end
                product = product * alpha;
            end
            sum_num = sum_num + product;
        end

        sum_num = sum_num + alpha_null( i );
        %sum_num is now the complete for row i

        %now to begin the denominator product
        for j = 1 : size( setsJ{ i }, 2 )
            jval = setsJ{ i }{ j };
            product_den = product_den * ( 1 + power( vecX( jval ) / matK(jval,i), matN(jval,i) ) );
        end

        new_vecX( i ) = sum_num / product_den;

    end
    
    if size( find( abs( vecX - new_vecX ) > 0.001 ), 1 ) == 0
        finished=true;
        vecX-new_vecX %display difference
    end

vecX = new_vecX;

end

fprintf('\nfinished\n');


