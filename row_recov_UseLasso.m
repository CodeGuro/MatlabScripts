function rec_row = row_recov_UseLasso( matDeltaX_selected_row, matDeltaX_selected, lambda )
    [ B, S ] = lasso( matDeltaX_selected', matDeltaX_selected_row' );
    rec_row = B(:, 1)';
    rec_row = matDeltaX_selected_row / matDeltaX_selected;
end
