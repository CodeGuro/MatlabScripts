% construct the vectors
n = 3; % size
matA = randi([1;10], n, n );
matA( logical( eye( n ) ) ) = 0;
vecb = randi( [0;10], length( matA ), 1 );
vecb_0 = vecb;
vecn = (rand( n, 1 ) - 0.5)/10;

% sanity check; find vector 'x' such that x = A*x + b
 vecx = ( eye( n ) - matA ) \ vecb;
 sanity_vecx = matA * vecx + vecb;

for itr = 0 : 10

    % find x
    vecx = ( eye( n ) - matA ) \ vecb;
    sanity_vecx = matA * vecx + vecb;

    matdeltaX = nan( n, n );

    % find the new states for all vecx(setdiff(1:n,i)) when vecx(i) is changed
    h = 1;
    for current = 1 : n
        selection = setdiff( 1:n, current );
        extracted_matA = matA( selection, current );
        % the arithmatic is performed here
        vecxp = (eye( n - 1 ) - matA( selection, selection ) )\( extracted_matA*( vecx(current)+h ) + vecb(selection) );
        vecxp = vecxp + itr * vecn(selection);
        % sanity check
        sanity_vecxp = matA( selection, selection ) * vecxp + extracted_matA*( vecx(current)+h ) + vecb( selection );

        %insert the x_i+h to vecxp to the i'th index
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
    vecx_noise{itr+1} = matA_rec * vecx + vecb;
    matA_recs_noise{itr+1}=matA_rec;
end