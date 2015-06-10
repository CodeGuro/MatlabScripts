function vecx_p = perturb_lin( n, matM, vecx, vecb, p_idx, p_amount, sigma )

    selection = setdiff( 1:n, p_idx );

    matM_p = matM;
    matM_p( p_idx, p_idx ) = 1 - (vecb(p_idx) / (vecx(p_idx)+p_amount));
    matM_p( p_idx, selection ) = 0;
    vecx_p = (eye(n) - matM_p) \ (vecb);
    vecx_p_sanity = matM_p*(vecx_p)+vecb;

    noise = randn(n, 1) * sigma;
    
    vecx_p = vecx_p + noise;

end

