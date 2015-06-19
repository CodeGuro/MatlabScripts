% generate the random n x n matrix
n = 20;
A_limiter = 0.8;
matA = rand(n,n) > A_limiter;
matA( logical( eye( n ) ) ) = 0;
matK = zeros( n, n );
matN = matA * 2;
vecX = zeros( n, 1 ); % the gene expressions
next_vecX = zeros( n, 1 ); % expression levels for next iteration
setsJ = {};
powerSetsJ = {};
alphas = {};
alpha_null = rand( n, 1 );

% generate setsJ, powerSetsJ, and matK
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
for i = 1 : size( powerSetsJ, 2 )
    for j = 1 : size( powerSetsJ{ i }, 2 )
        alphas{ i }{ j } = rand();
    end
end

%prepare for convergence
vecX = zeros( n, 1 );
iteration = 1;
maxIterations = 1000;
itDiff_threshold = 1E-8;
finished = false;


%start the convergence
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

% now that we've found a steady state, we can start the perturbations
perturb_amount = 1;
perturb_samples = 15;
steady_state_vecX = vecX;
steady_vecX = steady_state_vecX;
matdeltaX = NaN( n, n );
sigmas = 0:0.01:0.1;
vecMistakes_avg = nan(size(sigmas));
mistake_threshold = 1E-1;

for sigma_it = 1:size(sigmas,2)
    
    vecMistakes = nan( 1, perturb_samples );
    for sample_num = 1:perturb_samples
        
        for p = 1:n % perturbation index

            vecX = steady_vecX;
            vecX( p ) = steady_vecX( p ) + perturb_amount;
            iteration = 1;
            finished = false;

            while ~finished

                next_vecX = nonlinear_func( n, vecX, matK, matN, setsJ, powerSetsJ, alphas, alpha_null );
                next_vecX( p ) = steady_vecX( p ) + perturb_amount;

                if size( find( abs( vecX - next_vecX ) > itDiff_threshold ), 1 ) == 0 || iteration >= maxIterations
                    finished=true;
                    if iteration >= maxIterations
                        disp('warning! max iteration exceeded! aborting...');
                    end
                    vecDiff = vecX-next_vecX;

                    noise = randn(n, 1) * sigmas( sigma_it );
                    perturbed_i_steady_vecX = next_vecX;
                    perturbed_i_steady_vecX( p ) = steady_vecX( p ) + perturb_amount;
                    perturbed_i_steady_vecX = perturbed_i_steady_vecX + noise; % noise added here
                    vecdeltaX = perturbed_i_steady_vecX - steady_vecX;
                    matdeltaX( :, p ) = vecdeltaX;
                else
                    vecX = next_vecX;
                end
                iteration = iteration + 1;
            end
        end


        % We can now attempt to construct the linear matrix using the offsets
        % (deltamatX)
        lin_mat = nan( n, n );
        for current = 1 : n
            selection = setdiff( 1:n, current );
            lin_mat_cur_row = matdeltaX( current, selection ) / matdeltaX( selection, selection );
            lin_mat_cur_row = insert( lin_mat_cur_row, 0 , current );
            lin_mat( current, : ) = lin_mat_cur_row;
        end

        numMistakes = nnz( logical( matK ) - logical( abs(lin_mat) > mistake_threshold ) );
        vecMistakes( sample_num ) = numMistakes;
    end
    vecMistakes_avg( sigma_it ) = mean( vecMistakes );
end

plot( sigmas, vecMistakes_avg);
legend( 'avg mistakes' );
xlabel('sigma');
ylabel('recovery mistakes');
title(strcat(num2str(n),'X',num2str(n),'nonlinear matrix recovery errors with, ',num2str(perturb_samples),' samples per pertubation'));



