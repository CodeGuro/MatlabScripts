function [rec_row, rec_row_lambads] = row_recov_UseLasso( matDeltaX_selected_row, matDeltaX_selected, lambdas )
    [ B, S ] = lasso( matDeltaX_selected', matDeltaX_selected_row', 'lambda', lambdas );
    rec_row = B(:, 1)';
    rec_row_lambads = B';
end
