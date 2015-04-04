% %useful function: combnk(x,n) makes all combinations of n components from a
% vector 'x'

% % the set
% x = 1:4;
% for n = 1:4
% combnz = combnk(x,n);
% %or% combnz{nn} = combnk(x,nn);
% end

% generate the random 3x3 matrix
n = 5;
matA = randi( [0;1], n, n );
matA( logical( eye( n ) ) ) = 0;
matK = zeros( n, n );
matN = matA * 2; % since n_ij=1 where matA_ij=1
vecX = rand( [n, 1] ); % this will be refreshed with every iteration
new_vecX = zeros( n, 1 );

% fill the K(i,j) uniformly [0,1] where matA(i,j)==1
for i = 1 : n
    vecIdx = find( matA( i, : ) );
    for j = 1 : length( vecIdx )
        matK( i, vecIdx( j ) ) = rand();
    end
end

for repeat = 1 : 10000
    
    % for each row in matA, define J = { j | matA(i,j)==1 }
    for i = 1 : n
        vecJ = find( matA( i, : ) );
        sum_num = 0; % removing this comment causes problems % rand();
        product_den = 1; % product in denominator
        if size( vecJ ) < 1 
           % there are no nonzero elements for the current row
           new_vecX( i ) = 0;
           vecX( i ) = 0;
           continue; 
        end

        % find the power set of J
        powerSetJ = {};
        for it = 1 : size( vecJ,2 )
            powerSetJ{ it } = combnk( vecJ, it );
        end

        % begin evaluating the sum (numerator)
        % this loop resolves for the n'th element of power set of J. Ex: row [1,0,2]
        % yields {[1;3]},{[1,3]} which is equivalent to {{1},{3},{1,3}}
        % so we can treat rows as if they are subsets
        % [->1;3], [1;->3],   ->[1,3] evaluates to alpha_(it_idx)*Epsilon(1) +
        % alpha_(it_(idx+1))*Epsilon(3) + alpha_((it+1)_idx)*Epsilon(1)*Epsilon(3)
        for it = 1 : size( powerSetJ, 2 )
            vec_curElemPowerSetJ = powerSetJ{ it };
            % begin the sum (numerator) ex: [1,0,3] -> [1;3],
            for idx1 = 1 : size( vec_curElemPowerSetJ, 1 )
                product_num = 1; % product inside sum in numerator
                alpha = rand(); % ? for this subset of J
                for idx2 = 1 : size( vec_curElemPowerSetJ, 2 )
                  j = vec_curElemPowerSetJ( idx1, idx2 );
                  product_num = product_num * power( vecX(j) / matK(i,j), matN(i,j) );
                end
                sum_num = sum_num + alpha * product_num;
            end
        end

        %sum_num is now the complete for row i

        %now to begin the denominator product
        for it = 1:size( vecJ, 2 )
            j = vecJ( it );
            product_den = product_den * ( 1 + power( vecX(j) / matK(i,j), matN(i,j) ) );
        end

        new_vecX( i ) = sum_num / product_den;

    end
    
    vecX = new_vecX;
end

fprintf('\nfinished\n');


