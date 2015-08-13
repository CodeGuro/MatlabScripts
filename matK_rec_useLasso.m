function matK_recs = matK_rec_useLasso( n, matdeltaX, lambdas )
    matK_recs = nan( n, n, size(lambdas, 2) );
    
    for p_idx = 1 : n
        selection = setdiff( 1:n, p_idx );
        %matK_cur_rows = row_recov_UseLasso( matdeltaX( p_idx, selection ), matdeltaX( selection, selection ), lambdas );
        
        [ lasso_mat, nonzero_indices] = construct_perabolaMat( matdeltaX, n, p_idx );
        rowrecov = row_recov_UseLasso( matdeltaX( p_idx, selection ), lasso_mat, lambdas );
        matK_cur_rows = rowrecov( 1:length(lambdas), 1:(n-1) );
        matN_cur_rows = rowrecov( 1:length(lambdas), n:end );
        
        if( nnz(matN_cur_rows) > 0 )
            disp(matN_cur_rows);
        end
        
        Nmat = zeros( n, n );
        
        Nmat( nonzero_indices ) = matN_cur_rows( 1, :);
        
        insertion_element = 0;
        for z=1:size(matK_cur_rows,1)
            matK_recs( p_idx, :, z ) = insert(matK_cur_rows(z, :), insertion_element, p_idx);
        end
    end
    
end