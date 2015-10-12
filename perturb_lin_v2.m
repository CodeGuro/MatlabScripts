function vecx_p = perturb_lin_v2( n, steady_vecX, matK, vecb, p_idx, p_amount, sigma )

    selection = setdiff( 1:n, p_idx );
    extracted_matK = matK( selection, p_idx );
    % the arithmatic is performed here
    vecx_p = (eye( n - 1 ) - matK( selection, selection ) )\( extracted_matK*( steady_vecX(p_idx)+p_amount ) + vecb(selection) );
    % sanity check
    sanity_vecx_p = matK( selection, selection ) * vecx_p + extracted_matK*( steady_vecX(p_idx)+p_amount ) + vecb( selection );

    %insert the x_i+h to vecx_p to the i'th index
    noise = randn(n, 1) * sigma;
    noise( p_idx ) = 0;
    vecx_p = insert( vecx_p, steady_vecX( p_idx ) + p_amount, p_idx )';
    vecx_p = vecx_p + noise;

end

