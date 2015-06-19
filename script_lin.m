% construct the vectors
n = 10; % size
A_limiter = 0.1;
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

sigmas = 0:1E-3:0.1;
%pre allocate for speed
vecMistakes_avg = nan(size(sigmas));


% find the new states for all vecx(setdiff(1:n,i)) when vecx(i) is changed
for sigma_it = 1:size(sigmas,2)
    
    vecMistakes = nan( pert_samples, 1 );
    
    for sample_num = 1:pert_samples
        
        for p_idx = 1 : n
            vecx_p = perturb_lin_v2( n, matK, vecx, vecb, p_idx, pert_amount, sigmas(sigma_it) );
            vecdeltaX = vecx_p - vecx;

            %insert into the delta matrix
            matdeltaX( :, p_idx ) = vecdeltaX;
        end

        % we can now begin reconstruction of matK using the matdeltaX alone
        matK_rec = nan( n, n );
        for p_idx = 1 : n
            selection = setdiff( 1:n, p_idx );
            matK_cur_row = matdeltaX( p_idx, selection ) / matdeltaX( selection, selection );
            matK_cur_row = insert( matK_cur_row, 0 + sigmas( sigma_it ) * randn(), p_idx );
            matK_rec( p_idx, : ) = matK_cur_row;
        end
        
        numMistakes = nnz( logical( matK ) - logical( abs(matK_rec) > mistake_threshold ) );
        vecMistakes( sample_num ) = numMistakes;
        
    end

    vecMistakes_avg( sigma_it ) = mean( vecMistakes );

end


avg_smooth = smooth( sigmas, vecMistakes_avg, 0.25, 'rloess' );
plot( sigmas, vecMistakes_avg );
%plot( sigmas, vecMistakes_avg, sigmas, avg_smooth );
legend( 'avg mistakes' );
xlabel('sigma');
ylabel('recovery mistakes');
title(strcat(num2str(n),'X',num2str(n),' matrix recovery errors with, ',num2str(pert_samples),' samples per pertubation'));