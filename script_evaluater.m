load( 'ws_data_3.mat' );

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
        fig = figure();
        surf( perturbAmount_vec, mistakeThresh_real_vec, surface' );
        xlabel( 'perturbation amount' );
        ylabel( 'threshold' );
        zlabel( 'errors' );
        set(gca, 'XScale', 'log', 'YScale', 'log' );
        title( { 'perturbation X threshold ' ,[ 'slice at repetitions = ' num2str(numsamples_vec( it )) ] } );
        saveas( fig, [lin_sbisect_path '/rep_' num2str(it)] );
        close( fig );
    end
    
    % create the perturb x lasso surfaces
    lasso_sbisect_path = [ sbisect_path '/lasso' ];
    mkdir( lasso_sbisect_path );
    for it = 1:length(numsamples_vec)
        clear surface;
        surface( :, : ) = numMistakes_lasso( :, it, :, ms );
        fig = figure();
        surf( perturbAmount_vec, lambdas_vecreal, surface' );
        xlabel( 'perturbation amount' );
        ylabel( 'lambda' );
        zlabel( 'errors' );
        set(gca, 'XScale', 'log', 'YScale', 'log' );
        title( { 'perturbation X lambda ' ,[ 'slice at repetitions = ' num2str(numsamples_vec( it )) ] } );
        saveas( fig, [lasso_sbisect_path '/rep_' num2str(it)] );
        close( fig );
    end
    
    % create the perturb x qlasso surfaces
    qlasso_sbisect_path = [ sbisect_path '/Qlasso' ];
    mkdir( qlasso_sbisect_path );
    for it = 1:length(numsamples_vec)
        clear surface;
        surface( :, : ) = numMistakes_Qlasso( :, it, :, ms );
        fig = figure();
        surf( perturbAmount_vec, lambdas_vecreal, surface' );
        xlabel( 'perturbation amount' );
        ylabel( 'lambda' );
        zlabel( 'errors' );
        set(gca, 'XScale', 'log', 'YScale', 'log' );
        title( { 'perturbation X lambda : QLasso' ,[ 'slice at repetitions = ' num2str(numsamples_vec( it )) ] } );
        saveas( fig, [qlasso_sbisect_path '/rep_' num2str(it)] );
        close( fig );
    end
    
    % create the perturb x sample linear surface
    invsurf = min( numMistakes_inv, [], 3 );
    fig = figure();
    surf( perturbAmount_vec, numsamples_vec, invsurf( :, :, ms )' );
    xlabel( 'perturbation amount' );
    ylabel( 'samples' );
    zlabel( 'errors' );
    set(gca, 'XScale', 'log', 'YScale', 'log' );
    title( 'linear recovery method' );
    saveas( fig, [ surface_path '/surface_lin.fig' ] );
    close( fig );
    
    % create the perturb x sample lasso surface
    lassosurf = min( numMistakes_lasso, [], 3 );
    fig = figure();
    surf( perturbAmount_vec, numsamples_vec, lassosurf( :, :, ms )' );
    xlabel( 'perturbation amount' );
    ylabel( 'samples' );
    zlabel( 'errors' );
    set(gca, 'XScale', 'log', 'YScale', 'log' );
    title( 'linear recovery method' );
    saveas( fig, [ surface_path '/surface_lasso.fig' ] );
    close( fig );
    
    % create the perturb x sample Qlasso surface
    qlassosurf = min( numMistakes_Qlasso, [], 3 );
    fig = figure();
    surf( perturbAmount_vec, numsamples_vec, qlassosurf( :, :, ms )' );
    xlabel( 'perturbation amount' );
    ylabel( 'samples' );
    zlabel( 'errors' );
    set(gca, 'XScale', 'log', 'YScale', 'log' );
    title( 'linear recovery method' );
    saveas( fig, [ surface_path '/surface_Qlasso.fig' ] );
    close( fig );
    
    % create the repetitions comparison fig
    numsamp_mist(:,1) = min(min(numMistakes_inv(:,:,:,ms),[],3),[],1);
    numsamp_mist(:,2) = min(min(numMistakes_lasso(:,:,:,ms),[],3),[],1);
    numsamp_mist(:,3) = min(min(numMistakes_Qlasso(:,:,:,ms),[],3),[],1);
    fig = figure();
    plot( numsamples_vec', numsamp_mist );
    legend( 'linear', 'lasso', 'Qlasso' );
    xlabel( 'matrix repetitions' );
    ylabel( 'recovery errors' );
    set(gca, 'XScale', 'log' );
    title({'Error comparison for lambda, lasso, and Qlasso', ...
        'best-case across perturbations and lambda/threshold',['mat #' num2str(ms)]});
    saveas( fig, [ surface_path '/sample_comparison.fig' ] );
    close( fig );
    
    
end

% create repetitions comparison fig that averages all matrices
numsamp_mist(:,1) = min(min(mean(numMistakes_inv,4),[],3),[],1);
numsamp_mist(:,2) = min(min(mean(numMistakes_lasso,4),[],3),[],1);
numsamp_mist(:,3) = min(min(mean(numMistakes_Qlasso,4),[],3),[],1);
fig = figure();
plot( numsamples_vec', numsamp_mist );
legend( 'linear', 'lasso', 'Qlasso' );
xlabel( 'matrix repetitions' );
ylabel( 'recovery errors' );
set(gca, 'XScale', 'log' );
title({'Error comparison for lambda, lasso, and Qlasso', ...
    'best-case across perturbations and lambda/threshold','avg across all mats'});
saveas( fig, [fname_reps '/sample_comparison.fig' ] );
close( fig );



yvec_qlasso = nan( 1, length( perturbAmount_vec ) );
yvec_lasso = nan( 1, length( perturbAmount_vec ) );
yvec_inv = nan( 1, length( perturbAmount_vec ) );

for i=1:length( perturbAmount_vec )
    yvec_qlasso( i ) = min( min( numMistakes_Qlasso_avg(i,:,:) ) );
    yvec_lasso( i ) = min( min( numMistakes_lasso_avg(i,:,:) ) );
    yvec_inv( i ) = min( min( numMistakes_inv_avg(i,:,:) ) );
end

fig = figure();
plot( perturbAmount_vec, yvec_qlasso );
legend( 'errors across perturbations' );
xlabel( 'samples' ); ylabel( 'errors' );
title('min(lamdasxrepititions), avgMat, Quadratic lasso recovery smethod');
saveas( fig, [ fname_pert '/qlasso_avg' ] );
close( fig );

fig = figure();
plot( perturbAmount_vec, yvec_lasso );
legend( 'errors across perturbations' );
xlabel( 'samples' ); ylabel( 'errors' );
title('min(lamdasxrepititions), avgMat, lasso recovery smethod');
saveas( fig, [ fname_pert '/lasso_avg' ] );
close( fig );

fig = figure();
plot( perturbAmount_vec, yvec_inv );
legend( 'errors across perturbations' );
title('min(thresholdxrepititions), avgMat, linear recovery method');
xlabel( 'samples' ); ylabel( 'errors' );
saveas( fig, [ fname_pert '/lin_avg' ] );
close( fig );

% save bulk
qlasso_pert = [fname_pert '/qlasso'];
mkdir( qlasso_pert );
lasso_pert = [fname_pert '/lasso'];
mkdir( lasso_pert );
lin_pert = [fname_pert '/linear'];
mkdir( lin_pert );
for ms = 1:numMatSamplesExternal_int
    for i=1:length( perturbAmount_vec )
        yvec_qlasso( i ) = min( min( numMistakes_Qlasso(i,:,:,ms) ) );
        yvec_lasso( i ) = min( min( numMistakes_lasso(i,:,:,ms) ) );
        yvec_inv( i ) = min( min( numMistakes_inv(i,:,:,ms) ) );
    end
    
    fig = figure();
    plot( perturbAmount_vec, yvec_qlasso );
    legend( 'errors across perturbations' );
    xlabel( 'perturbation' ); ylabel( 'errors' );
    title({'min(lamdasxrepititions), quadratic lasso recovery smethod', ['mat #' num2str(ms)]});
    saveas( fig, [ qlasso_pert '/lasso_mat_' num2str(ms) ] );
    close( fig );
    
    fig = figure();
    plot( perturbAmount_vec, yvec_lasso );
    legend( 'errors across perturbations' );
    xlabel( 'perturbation' ); ylabel( 'errors' );
    title({'min(lamdasxrepititions), lasso recovery smethod', ['mat #' num2str(ms)]});
    saveas( fig, [ lasso_pert '/lasso_mat_' num2str(ms) ] );
    close( fig );

    fig = figure();
    plot( perturbAmount_vec, yvec_inv );
    legend( 'errors across perturbations' );
    title({'min(thresholdxrepititions), linear recovery method', ['mat #' num2str(ms)]});
    xlabel( 'perturbation' ); ylabel( 'errors' );
    saveas( fig, [ lin_pert '/inv_mat_' num2str(ms) ] );
    close( fig );
end
