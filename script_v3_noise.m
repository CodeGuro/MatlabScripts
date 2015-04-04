% construct the vectors
n = 3; % size
matA = randi([1;10], n, n );
matA( logical( eye( n ) ) ) = 0;
vecb = randi( [0;10], length( matA ), 1 );
vecb_0 = vecb;
itmax = 100;


% sanity check; find vector 'x' such that x = A*x + b
 vecx = ( eye( n ) - matA ) \ vecb;
 sanity_vecx = matA * vecx + vecb;
 

for itr = 0 : itmax

    matNoise = randn( n, n ); % will be scaled later
    sigma = itr/itmax * 1;
    sigmavec(itr+1)=sigma;
    % find x
    vecx = ( eye( n ) - matA ) \ vecb;
    sanityCheck = matA * vecx + vecb;
    
    % x' = x + noise
    vecx = vecx + sigma * randn(n,1);

    matdeltaX = nan( n, n );

    % find the new states after perturbing vecx(i)
    h = 1;
    for current = 1 : n
        selection = setdiff( 1:n, current );
        extracted_matA = matA( selection, current );
        % the arithmatic is performed here
        vecxp = (eye( n - 1 ) - matA( selection, selection ) )\( extracted_matA*( vecx(current)+h ) + vecb(selection) );

        %insert the x_i' to vecxp at the i'th index
        vecxp = insert( vecxp, vecx( current ) + h, current )';
        vecdeltaX = vecxp - vecx;

        %insert into the delta matrix
        matdeltaX( current, : ) = vecdeltaX;
    end

    % we can now begin reconstruction of matA using the matdeltaX alone
    matA_rec = nan( n, n );
    for current = 1 : n
        selection = setdiff( 1:n, current );
        matA_cur_row = matdeltaX( selection, selection ) \ matdeltaX( selection, current );
        matA_cur_row = insert( matA_cur_row, 0, current );
        matA_rec( current, : ) = matA_cur_row;
    end
    
    for i = 1 : n
        for j = 1 : n
            matA_rec_noise_plots((i-1)*n + j, itr+1) = matA_rec(i,j);
        end
        vecx_noise_plots(:,itr+1)= matA_rec*vecx+vecb;
    end
end


for i = 1 : n*n
	plot(sigmavec,matA_rec_noise_plots(i,:))
    xlabel('sigma');
    mat_i = floor((i-1)/n);
    mat_j = i-(n*mat_i);
    ylabel(strcat('matA_{',int2str(mat_i+1),int2str(mat_j),'}'));
   % axis([ 0 sigma -30 30 ] );
    pause(2)
end


for i=1:n
    plot(sigmavec,vecx_noise_plots(i,:));
    xlabel('sigma');
    ylabel(strcat('x_',int2str(i)));
   % axis([ 0 sigma -30 30 ] )
    pause(2);
end