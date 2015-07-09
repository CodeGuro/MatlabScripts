% generate the random n x n matrix
n = 5;
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
itDiff_threshold = 1E-4;
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
perturb_samples = 20;
steady_vecX = vecX;
matdeltaX = NaN( n, n );
sigmas = 0:1E-2:0.1;
vecMistakes_avg = nan(size(sigmas));
vecstDevs = nan(size(sigmas));
mistake_threshold = 1E-1;
lambdas = [2E-3 2E-2 2E-1 2E0 2E1 ];
vecMistakes_avg_lasso = nan( length(lambdas), length(sigmas) );
vecstDevs_lasso = nan( length(lambdas), length(sigmas) );

for sigma_it = 1:length(sigmas)
    disp(['sampling sigma: ' num2str(sigma_it) ' out of ' num2str(length(sigmas))]);
    vecMistakes = nan( 1, perturb_samples );
    vecMistakes_lasso = nan( perturb_samples, length(lambdas) );
    
    for sample_num = 1:perturb_samples
        
        disp(['iterating... sample ' num2str(sample_num) ' out of ' num2str(perturb_samples)]);
        
        for p = 1:n % perturbation index

            disp(['perturbation index: ' num2str(p) ' out of ' num2str(n)]);
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
                    noise(p) = 0;
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
        matK_rec = matK_rec_useInv( n, matdeltaX );
        numMistakes = nnz( logical( matK ) - logical( abs(matK_rec) > mistake_threshold ) );
        vecMistakes( sample_num ) = numMistakes;
        
        if length(lambdas) > 0
            matK_recs_lasso = matK_rec_useLasso( n, matdeltaX, lambdas );
            for z=1:length(lambdas)
                numMistakes_lasso = nnz( logical( matA ) - logical( abs( matK_recs_lasso(:,:,z) ) > mistake_threshold ) );
                vecMistakes_lasso( sample_num, z ) = numMistakes_lasso;
            end
        end
        
    end
    vecMistakes_avg( sigma_it ) = mean( vecMistakes );
    vecstDevs( sigma_it ) = std( vecMistakes );
    
    if length(lambdas) > 0
        vecMistakes_avg_lasso( :, sigma_it ) = mean( vecMistakes_lasso );
        vecstDevs_lasso( :, sigma_it ) = std( vecMistakes_lasso );
    end
    
end

% output
plot_sigmas = (ones(length(lambdas),1)*sigmas)';
plot_mistakes = vecMistakes_avg_lasso';
plot_devs = vecstDevs_lasso';
legends = cell( 1, length(lambdas) );
for lambda_it = 1:length(lambdas)
    legends{ lambda_it } = ['lasso: lambda=' num2str( lambdas(lambda_it) ) ];
end

plot_sigmas( :, end+1 ) = sigmas;
plot_mistakes(:, end+1 ) = vecMistakes_avg;
plot_devs(:, end+1 ) = vecstDevs;
legends{length(lambdas) + 1} = 'inv';

errorbar( plot_sigmas, plot_mistakes, plot_devs, ':o' );
legend( 'avg mistakes' );
xlabel('sigma');
ylabel('recovery mistakes');
legend( legends );
title( {[num2str(n) 'X' num2str(n) ' matrix of a linear regulatory network, with '...
    num2str(length(sigmas)) ' samples, with '...
    num2str(perturb_samples) ' samples per perturbation'],...
    ['perturb amount=' num2str(perturb_amount) ', '...
    'A limiter=' num2str(A_limiter) ', ' 'nnz(matA)=' num2str(nnz(matA))]} );
%axis( [ sigmas(1) sigmas(end) 0 10] );
axis('auto');
