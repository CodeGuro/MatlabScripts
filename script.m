% construct the vectors
n = 10; % size
matA = randi([1;10], n, n );
matA( logical( eye( n ) ) ) = 0;
b_min = 1E-5;
b_max = 10;
vecb = b_min + (b_max - b_min).*rand( [n,1] );
vecMistakes = [];
sigmas = 0:1E-3:0.1;

% find vector 'x' such that x = A*x + b
vecx = ( eye( n ) - matA ) \ vecb;
sanity_vecx = matA * vecx + vecb;
matdeltaX = nan( n, n );

% find the new states for all vecx(setdiff(1:n,i)) when vecx(i) is changed
for sigma_it = 1:size(sigmas,2)
    h = 1;
    
    for current = 1 : n
        vecx_p = perturb_lin( n, matA, vecx, vecb, current, h, sigmas(sigma_it) );
        vecdeltaX = vecx_p - vecx;

        %insert into the delta matrix
        matdeltaX( current, : ) = vecdeltaX;
    end

    % we can now begin reconstruction of matA using the matdeltaX alone
    matA_rec = nan( n, n );
    for current = 1 : n
        selection = setdiff( 1:n, current );
        matA_cur_row = matdeltaX( selection, selection ) \ matdeltaX( selection, current );
        matA_cur_row = insert( matA_cur_row, 0 + sigmas( sigma_it ) * randn(), current );
        matA_rec( current, : ) = matA_cur_row;
    end

    numMistakes = abs( nnz( matA ) - nnz( matA_rec > 1E-3) );
    vecMistakes( sigma_it, 1 ) = numMistakes;
end

ys = smooth( sigmas, vecMistakes, 0.25, 'rloess' );
plot( sigmas, vecMistakes, sigmas, ys );
legend( 'raw sample data', 'smoothed samle data' );
xlabel('sigma');
ylabel('recov_mistakes');