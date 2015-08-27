function matK_recs = matK_rec_useLasso( n, matdeltaX, num_samples, lambdas, use_lasso_Nmat )
    matK_recs = nan( n, n, size(lambdas, 2) );
    
    for p_idx = 1 : n
        selection_rows = setdiff( 1:n, p_idx );
        selection_cols = setdiff( 1:n*num_samples, p_idx + (0:(num_samples-1))*n );
        
        if use_lasso_Nmat
            
            lasso_mat = nan( n*(n+1)/2 - 1, (n-1)*num_samples );
            for num_sample = 1:num_samples
                lasso_mat( :, ((num_sample-1)*(n-1)+1):num_sample*(n-1) )...
                    = construct_perabolaMat( matdeltaX( :, ((num_sample-1)*n+1):num_sample*n ), n, p_idx );
            end
            
            [ placeholder, nonzero_indices ] = construct_perabolaMat( matdeltaX, n, p_idx );
            rowrecov = row_recov_UseLasso( matdeltaX( p_idx, selection_cols ), lasso_mat, lambdas );
            matK_cur_rows = rowrecov( 1:length(lambdas), 1:(n-1) );
            matN_cur_rows = rowrecov( 1:length(lambdas), n:end );
        else
            matK_cur_rows = row_recov_UseLasso( matdeltaX( p_idx, selection_cols ), matdeltaX( selection_rows, selection_cols ), lambdas );
        end
        
%         if( nnz(matN_cur_rows) > 0 )
%             disp(matN_cur_rows);
%         end
        
        insertion_element = 0;
        for z=1:size(matK_cur_rows,1)
            matK_recs( p_idx, :, z ) = insert(matK_cur_rows(z, :), insertion_element, p_idx);
        end
        
        if use_lasso_Nmat
            for lambda_idx = 1:length(lambdas)
                Nmat_gene_pidx = zeros( n, n );
                Nmat_gene_pidx( nonzero_indices ) = matN_cur_rows( lambda_idx, : );
                
         %       [ nnz_I, nnz_J ] = find( Nmat_gene_pidx );
         %       matK_recs( p_idx, nnz_I, lambda_idx ) = Inf;
         %       matK_recs( p_idx, nnz_J, lambda_idx ) = Inf;
                
                for i=1:n
                    for j=1:n
                        if Nmat_gene_pidx( i, j ) * Nmat_gene_pidx( i, j ) > 1E-2
                            matK_recs( p_idx, i, lambda_idx ) = Inf;
                            matK_recs( p_idx, j, lambda_idx ) = Inf;
                        end
                    end
                end
            end
        end
        
    end
    
end