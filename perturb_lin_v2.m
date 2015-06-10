function vecx_p = perturb_lin_v2( n, matM, vecx, vecb, p_idx, p_amount, sigma )

    selection = setdiff( 1:n, p_idx );
    extracted_matM = matM( selection, p_idx );
    % the arithmatic is performed here
    vecx_p = (eye( n - 1 ) - matM( selection, selection ) )\( extracted_matM*( vecx(p_idx)+p_amount ) + vecb(selection) );
    % sanity check
    sanity_vecx_p = matM( selection, selection ) * vecx_p + extracted_matM*( vecx(p_idx)+p_amount ) + vecb( selection );

    %insert the x_i+h to vecx_p to the i'th index
    noise = randn(n, 1) *sigma;
    vecx_p = insert( vecx_p, vecx( p_idx ) + p_amount, p_idx )';
    vecx_p = vecx_p + noise;

end

