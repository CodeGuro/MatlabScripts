function matK_recs = matK_rec_useLasso( n, matdeltaX, lambdas )

    matK_rec = nan( n, n );
    matK_recs = nan( n, n, size(lambdas, 2) );
    
    for p_idx = 1 : n
        selection = setdiff( 1:n, p_idx );
        [ matK_cur_row , matK_cur_rows ] = row_recov_UseLasso( matdeltaX( p_idx, selection ), matdeltaX( selection, selection ), lambdas );
        
        insertion_element = 0;
        matK_rec( p_idx, : ) = insert( matK_cur_row, insertion_element, p_idx );
        
        for z=1:size(matK_cur_rows,1)
            matK_recs( p_idx, :, z ) = insert(matK_cur_rows(z, :), insertion_element, p_idx);
        end
    end
    
end