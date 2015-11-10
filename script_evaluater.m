load( 'ws_data_2.mat' );

fname_pert = 'perturbations';
fname_reps = 'repetitions';

mkdir( fname_pert );
mkdir( fname_reps );

% save bulk
for ms = 1:numMatSamplesExternal_int
    
    mat_path = [ fname_reps '/mat' num2str(ms) ];
    mkdir( mat_path );
    
    surface_path = [mat_path '/surface'];
    mkdir( surface_path );

    
    % create the perturb x thresh surfaces
    sbisect_path = [ surface_path '/surface_bisect' ];
    mkdir( sbisect_path );
    
    lin_sbisect_path = [ sbisect_path '/linear' ];
    mkdir( lin_sbisect_path );
    for it = 1:length(numsamples_vec)
        clear surface;
        surface( :, : ) = numMistakes_inv( :, it, :, ms );
        fig = surf( perturbAmount_vec, mistakeThresh_real_vec, surface' );
        xlabel( 'perturbation amount' );
        ylabel( 'threshold' );
        zlabel( 'mistakes' );
        set(gca, 'XScale', 'log', 'YScale', 'log' );
        title( { 'perturbation X threshold ' ,[ 'slice at repetitions = ' num2str(numsamples_vec( it )) ] } );
        saveas( fig, [lin_sbisect_path '/rep_' num2str(it)] );
    end
    
    % create the perturb x lasso surfaces
    lasso_sbisect_path = [sbisect_path '/lasso' ];
    mkdir( lasso_sbisect_path );
    for it = 1:length(lambdas_vecreal)
        clear surface;
        surface( :, : ) = numMistakes_lasso( :, it, :, ms );
        fig = surf( perturbAmount_vec, lambdas_vecreal, surface' );
        xlabel( 'perturbation amount' );
        ylabel( 'lambda' );
        zlabel( 'mistakes' );
        set(gca, 'XScale', 'log', 'YScale', 'log' );
        title( { 'perturbation X lambda ' ,[ 'slice at repetitions = ' num2str(numsamples_vec( it )) ] } );
        saveas( fig, [lasso_sbisect_path '/rep_' num2str(it)] );
    end
    
    % create the perturb x sample linear surface
    invsurf = min( numMistakes_inv, [], 3 );
    fig = surf( perturbAmount_vec, numsamples_vec, invsurf( :, :, ms )' );
    xlabel( 'perturbation amount' );
    ylabel( 'samples' );
    zlabel( 'mistakes' );
    set(gca, 'XScale', 'log', 'YScale', 'log' );
    title( 'linear recovery method' );
    saveas( fig, [ surface_path '/surface_lin.fig' ] );
    
    % create the perturb x sample lasso surface
    lassosurf = min( numMistakes_lasso, [], 3 );
    fig = surf( perturbAmount_vec, numsamples_vec, lassosurf( :, :, ms )' );
    xlabel( 'perturbation amount' );
    ylabel( 'samples' );
    zlabel( 'mistakes' );
    set(gca, 'XScale', 'log', 'YScale', 'log' );
    title( 'linear recovery method' );
    saveas( fig, [ surface_path '/surface_lasso.fig' ] );
    
    
end

yvec_lasso = nan( 1, length( perturbAmount_vec ) );
yvec_inv = nan( 1, length( perturbAmount_vec ) );

for i=1:length( perturbAmount_vec )
    yvec_lasso( i ) = min( min( numMistakes_lasso_avg(i,:,:) ) );
    yvec_inv( i ) = min( min( numMistakes_inv_avg(i,:,:) ) );
end

fig = plot( perturbAmount_vec, yvec_lasso );
legend( 'mistakes across perturbations' );
xlabel( 'samples' ); ylabel( 'mistakes' );
title('min(lamdasxrepititions), avgMat, lasso recovery smethod');
saveas( fig, [ fname_pert '/lasso_avg' ] );

fig = plot( perturbAmount_vec, yvec_inv );
legend( 'mistakes across perturbations' );
title('min(thresholdxrepititions), avgMat, linear recovery method');
xlabel( 'samples' ); ylabel( 'mistakes' );
saveas( fig, [ fname_pert '/lin_avg' ] );

% save bulk
lasso_pert = [fname_pert '/lasso'];
mkdir( lasso_pert );
lin_pert = [fname_pert '/linear'];
mkdir( lin_pert );
for ms = 1:numMatSamplesExternal_int
    for i=1:length( perturbAmount_vec )
        yvec_lasso( i ) = min( min( numMistakes_lasso(i,:,:,ms) ) );
        yvec_inv( i ) = min( min( numMistakes_inv(i,:,:,ms) ) );
    end
    
    fig = plot( perturbAmount_vec, yvec_lasso );
    legend( 'mistakes across perturbations' );
    xlabel( 'perturbation' ); ylabel( 'mistakes' );
    title({'min(lamdasxrepititions), lasso recovery smethod', ['mat #' num2str(ms)]});
    saveas( fig, [ lasso_pert '/lasso_mat_' num2str(ms) ] );

    fig = plot( perturbAmount_vec, yvec_inv );
    legend( 'mistakes across perturbations' );
    title({'min(thresholdxrepititions), linear recovery method', ['mat #' num2str(ms)]});
    xlabel( 'perturbation' ); ylabel( 'mistakes' );
    saveas( fig, [ lin_pert '/inv_mat_' num2str(ms) ] );
end
