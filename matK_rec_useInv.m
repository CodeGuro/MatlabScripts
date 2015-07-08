function matK_rec = matK_rec_useInv( n, matdeltaX )

    matK_rec = nan( n, n );
    
    for p_idx = 1 : n
        selection = setdiff( 1:n, p_idx );
        matK_cur_row = row_recov_UseInv( matdeltaX( p_idx, selection ), matdeltaX( selection, selection ) );
        matK_cur_row = insert( matK_cur_row, 0, p_idx );
        matK_rec( p_idx, : ) = matK_cur_row;
    end
    
end