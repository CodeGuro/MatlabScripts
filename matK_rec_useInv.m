function matK_rec = matK_rec_useInv( n, matdeltaX, num_samples )

    matK_rec = nan( n, n );
    
    for p_idx = 1 : n
        selection_rows = setdiff( 1:n, p_idx );
        selection_cols = setdiff( 1:n*num_samples, p_idx + (0:(num_samples-1))*n );
        matK_cur_row = row_recov_UseInv( matdeltaX( p_idx, selection_cols ), matdeltaX( selection_rows, selection_cols ) );
        matK_cur_row = insert( matK_cur_row, 0, p_idx );
        matK_rec( p_idx, : ) = matK_cur_row;
    end
    
end