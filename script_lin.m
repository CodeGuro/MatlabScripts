% construct the vectors
n = 10; % size
A_limiter = 0.8;
mistake_threshold = 1E-1;
lb = mistake_threshold;
ub = 1;
matA = rand(n,n) > 0.8;
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
pert_samples = 50;

sigmas = 0:1E-3:0.1;
%pre allocate for speed
vecMistakes_avg = nan(size(sigmas));
vecDevMin_avg = nan(size(sigmas));
vecDevMax_avg = nan(size(sigmas));


% find the new states for all vecx(setdiff(1:n,i)) when vecx(i) is changed
for sigma_it = 1:size(sigmas,2)
    
    vecMistakes = nan( pert_samples, 1 );
    min_devs = nan( pert_samples, 1 );
    max_devs = nan( pert_samples, 1 );
    
    for sample_num = 1:pert_samples
        
        for p_idx = 1 : n
            vecx_p = perturb_lin_v2( n, matK, vecx, vecb, p_idx, pert_amount, sigmas(sigma_it) );
            vecdeltaX = vecx_p - vecx;

            %insert into the delta matrix
            matdeltaX( p_idx, : ) = vecdeltaX;
        end

        % we can now begin reconstruction of matK using the matdeltaX alone
        matK_rec = nan( n, n );
        for p_idx = 1 : n
            selection = setdiff( 1:n, p_idx );
            matK_cur_row = matdeltaX( selection, selection ) \ matdeltaX( selection, p_idx );
            matK_cur_row = insert( matK_cur_row, 0 + sigmas( sigma_it ) * randn(), p_idx );
            matK_rec( p_idx, : ) = matK_cur_row;
        end
        
        numMistakes = nnz( logical( matK ) - logical( abs(matK_rec) > mistake_threshold ) );
        vecMistakes( sample_num ) = numMistakes;
        min_devs( sample_num ) = min(vecMistakes);
        max_devs( sample_num ) = max(vecMistakes);
        
    end

    vecMistakes_avg( sigma_it ) = mean( vecMistakes );
    vecDevMin_avg( sigma_it ) = mean( min_devs );
    vecDevMax_avg( sigma_it ) = mean( max_devs );
    
end


avg_smooth = smooth( sigmas, vecMistakes_avg, 0.25, 'rloess' );
devMin_smooth = smooth( sigmas, vecDevMin_avg, 0.25, 'rloess' );
devMax_smooth = smooth( sigmas, vecDevMax_avg, 0.25, 'rloess' );
plot( sigmas, vecMistakes_avg, sigmas, avg_smooth, sigmas, vecDevMin_avg, sigmas, devMin_smooth, sigmas, vecDevMax_avg, sigmas, devMax_smooth );
%plot( sigmas, vecMistakes_avg, sigmas, avg_smooth );
legend( 'avg mistakes (raw)', 'avg mistakes (smooth)', 'avg min deviation (raw)', 'avg min deviation (smooth)', 'avg max deviation (raw)', 'avg max deviation (smooth)' );
xlabel('sigma');
ylabel('recovery mistakes');
title(strcat(num2str(n),'X',num2str(n),' matrix recovery errors with, ',num2str(pert_samples),' samples per pertubation'));