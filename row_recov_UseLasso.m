function rec_row = row_recov_UseLasso( matDeltaX_selected_row, matDeltaX_selected, lambda )
    rec_row_orig = row_recov_UseInv( matDeltaX_selected_row, matDeltaX_selected );
    [ B, S ] = lasso( matDeltaX_selected', matDeltaX_selected_row', 'lambda', lambda );
    rec_row = B(:, 1)';
end
