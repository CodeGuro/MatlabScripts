% construct the vectors
n = 10; % size
A_limiter = 0.3;
mistake_threshold = 1E-1;
lb = mistake_threshold;
ub = 1;
matA = rand(n,n) > A_limiter;
matA( logical( eye( n ) ) ) = 0;
matK = zeros(n,n); %declare so dimensions match
matK(matA) = lb + (ub-lb).*rand(nnz(matA),1); % K_ij =/= 0 --> A_ij =/= 0
b_min = 1E-5;
b_max = ub;
vecb = b_min + (b_max - b_min).*rand( [n,1] );

% find vector 'x' such that x = A*x + b
vecx = ( eye( n ) - matK ) \ vecb;
sanity_vecx = matK * vecx + vecb;

%prepare for pertubations
matdeltaX = nan( n, n );
pert_amount = 1; % this should be scaled up for best results on average
pert_samples = 20;
perturb_func = @perturb_lin_v2;

sigmas = 0:0.01:0.3;
lambdas = [2E-1 2E-2 2E-3 2E-4 2E-5];
%pre allocate for speed
vecMistakes_avg = nan(size(sigmas));
vecstDevs = nan(size(sigmas));

vecMistakes_avg_lasso = nan( size(lambdas,2), size(sigmas,2) );
vecstDevs_lasso = nan( size(lambdas,2), size(sigmas,2) );


% find the new states for all vecx(setdiff(1:n,i)) when gene i is changed
for sigma_it = 1:size(sigmas,2)
    disp(['sampling sigma: ' num2str(sigma_it) ' out of ' num2str(size(sigmas,2))]);
    
    vecMistakes = nan( pert_samples, 1 );
    vecMistakes_lasso = nan( pert_samples, size( lambdas, 2 ) );
    
    for sample_num = 1:pert_samples
        
        disp(['iterating... sample ' num2str(sample_num) ' out of ' num2str(pert_samples)]);
        
        for p_idx = 1 : n
            vecx_p = perturb_func( n, matK, vecx, vecb, p_idx, pert_amount, sigmas(sigma_it) );
            vecdeltaX = vecx_p - vecx;

            %insert into the delta matrix
            matdeltaX( :, p_idx ) = vecdeltaX;
        end

        % we can now begin reconstruction of matK using the matdeltaX alone
        matK_rec = matK_rec_useInv( n, matdeltaX, sigmas, sigma_it );
        numMistakes = nnz( logical( matK ) - logical( abs(matK_rec) > mistake_threshold ) );
        vecMistakes( sample_num ) = numMistakes;
        
       	matK_recs_lasso = matK_rec_useLasso( n, matdeltaX, sigmas, sigma_it, lambdas );
        for z=1:size(lambdas,2)
            numMistakes_lasso = nnz( logical( matK_recs_lasso(:,:,z) ) - logical( abs(matK_rec) > mistake_threshold ) );
            vecMistakes_lasso( sample_num, z ) = numMistakes_lasso;
        end
        
    end

    vecMistakes_avg( sigma_it ) = mean( vecMistakes );
    vecstDevs( sigma_it ) = std( vecMistakes );
    
    vecMistakes_avg_lasso( :, sigma_it ) = mean( vecMistakes_lasso );
    vecstDevs_lasso( :, sigma_it ) = std( vecMistakes_lasso );

end

plot_sigmas = (ones(size(lambdas,2),1)*sigmas)';
plot_mistakes = vecMistakes_avg_lasso';
plot_devs = vecstDevs_lasso';
for lambda_it = 1:size(lambdas,2)
    legends{lambda_it} = num2str( lambdas(lambda_it) );
end

plot_sigmas( :, end+1 ) = sigmas;
plot_mistakes(:, end+1 ) = vecMistakes_avg;
plot_devs(:, end+1 ) = vecstDevs;
legends{size(lambdas,2) + 1} = 'inv';

errorbar( plot_sigmas, plot_mistakes, plot_devs, ':o' );
legend( 'avg mistakes' );
xlabel('sigma');
ylabel('recovery mistakes');
legend( legends );
title( {[num2str(n) 'X' num2str(n) ' matrix of a linear regulatory network, with '...
    num2str(length(sigmas)) ' samples, with '...
    num2str(pert_samples) ' samples per perturbation'],...
    ['perturb amount=' num2str(pert_amount) ', '...
    'A limiter=' num2str(A_limiter)]} );
%axis( [ sigmas(1) sigmas(end) 0 10] );
axis('auto');