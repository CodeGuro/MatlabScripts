function result = construct_perabolaMat( matdeltaX, n, i )
    
    % N_ijk -> gene j and k effects on gene i
    
    % Construct the unprocessed matrix
    rowCount = n*(n+1)/2;
    columnCount = n;
    peraMat = nan( rowCount, columnCount );

    for segment = 1:n
        i_para = rowCount - (n-(segment-1))*(n-(segment-2))/2 + 1;
        for row = segment:n
            for column = 1:n
                j_para = column;
                peraMat( i_para, j_para ) = matdeltaX( segment, column ) * matdeltaX( row, column );
            end
            i_para = i_para + 1;
        end
    end
    
    % knock out column i (removes redundant terms), N_iji=0
    column_selection = setdiff( 1:n, i );
    
    % knock out n rows (removes redundant terms), N_iij=0
    row_selection = 1:(n*(n+1)/2);
    n_selections = reshape( 1:(n*n), [n n] )';
    for k=1:n
        n_selections(k,:) = n_selections(k,:) - (k-1)*k/2;
        n_selections( k, 1:(k-1) ) = 0;
    end
    
    n_selections( i, : ) = 0;
    n_selections( :, i ) = 0;
    
    row_selection = setdiff( reshape( n_selections', [ 1 n*n ] ), 0 );
    
    result = [ matdeltaX( column_selection, column_selection );...
        peraMat( row_selection, column_selection )];
    
end

