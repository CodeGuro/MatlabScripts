load( 'ws_data_2.mat' );

fname_pert = 'perturbations';
fname_reps = 'repetitions';

mkdir( fname_pert );
mkdir( fname_reps );

yvec_lasso = nan( 1, length( numsamples_vec ) );
yvec_inv = nan( 1, length( numsamples_vec ) );

for j=1:length( numsamples_vec )
    yvec_lasso( j ) = min( min( numMistakes_lasso_avg(:,j,:) ) );
    yvec_inv( j ) = min( min( numMistakes_inv_avg(:,j,:) ) );
end

fig = plot( numsamples_vec, yvec_lasso );
legend( 'mistakes across repetitions' );
xlabel( 'samples' ); ylabel( 'mistakes' );
title('min(lambdasXperturb), avgMat, lasso recovery smethod');
saveas( fig, [ fname_reps '/lasso_avg' ] );

fig = plot( numsamples_vec, yvec_inv );
legend( 'mistakes across repetitions' );
title('min(thresholdxperturb), avgMat, linear recovery method');
xlabel( 'samples' ); ylabel( 'mistakes' );
saveas( fig, [fname_reps '/inv_avg' ] );

% save bulk
for ms = 1:numMatSamplesExternal_int
    for j=1:length( numsamples_vec )
        yvec_lasso( j ) = min( min( numMistakes_lasso(:,j,:,ms) ) );
        yvec_inv( j ) = min( min( numMistakes_inv(:,j,:,ms) ) );
    end
    
    fig = plot( numsamples_vec, yvec_lasso );
    legend( 'mistakes across repetitions' );
    xlabel( 'samples' ); ylabel( 'mistakes' );
    title('min(lambdasXperturb), avgMat, lasso recovery smethod');
    saveas( fig, [ fname_reps '/lasso_sample_' num2str(ms) ] );

    fig = plot( numsamples_vec, yvec_inv );
    legend( 'mistakes across repetitions' );
    title('min(thresholdxperturb), avgMat, linear recovery method');
    xlabel( 'samples' ); ylabel( 'mistakes' );
    saveas( fig, [ fname_reps '/inv_sample_' num2str(ms) ] );
    
    mat_path = [ fname_reps '/mat' num2str(ms) ];
    mkdir( mat_path );
    bisect_path = [ mat_path '/bisect' ];
    mkdir( bisect_path );
  %{  
    for perturbidx = 1:length( perturbAmount_vec )
        perturb_path = [ bisect_path '/perturb_' num2str(perturbidx)];
        mkdir(perturb_path);
       
        thresh_path = [ perturb_path '/thresh' ];
        mkdir(thresh_path);
        for threshidx = 1:length(mistakeThresh_real_vec)
            yvec_inv_thresh = numMistakes_inv( perturbidx, :, threshidx, ms );
            fig = plot( numsamples_vec, yvec_inv_thresh );
            title( [ 'threshold = ' num2str( mistakeThresh_real_vec(threshidx) )...
                ' at perturb=' num2str(perturbAmount_vec(perturbidx)) ] );
            xlabel( 'repetitions' );
            ylabel( 'recov mistakes' );
            fname = [ thresh_path '/thresh_' num2str( threshidx ) ];
            saveas( fig, fname );
            disp(['saving: ' fname] );
        end
        
        lambda_path = [ perturb_path '/lambda' ];
        mkdir(lambda_path);
        for lambdaidx = 1:length(lambdas_vecreal)
            yvec_lasso_lambda = numMistakes_lasso( perturbidx, :, lambdaidx, ms );
            fig = plot( numsamples_vec, yvec_lasso_lambda );
            title( [ 'lambda = ' num2str( lambdas_vecreal(lambdaidx) )...
                ' at perturb=' num2str(perturbAmount_vec(perturbidx)) ] );
            xlabel( 'repetitions' );
            ylabel( 'recov mistakes' );
            fname = [ lambda_path '/lambda_' num2str( lambdaidx ) ];
            saveas( fig, fname );
            disp( ['saving: ' fname ] );
        end        
    end
    %}
    
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
        title( { 'perturbation X lambda ' ,[ 'slice at repetitions = ' num2str(numsamples_vec( it )) ] } );
        saveas( fig, [lasso_sbisect_path '/rep_' num2str(it)] );
    end
    
    % create the perturb x sample linear surface
    invsurf = min( numMistakes_inv, [], 3 );
    fig = surf( perturbAmount_vec, numsamples_vec, invsurf( :, :, ms )' );
    xlabel( 'perturbation amount' );
    ylabel( 'samples' );
    zlabel( 'mistakes' );
    title( 'linear recovery method' );
    saveas( fig, [ surface_path '/surface_lin.fig' ] );
    % create the perturb x sample lasso surface
    lassosurf = min( numMistakes_lasso, [], 3 );
    fig = surf( perturbAmount_vec, numsamples_vec, lassosurf( :, :, ms )' );
    xlabel( 'perturbation amount' );
    ylabel( 'samples' );
    zlabel( 'mistakes' );
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
saveas( fig, [ fname_pert '/inv_avg' ] );

% save bulk
for ms = 1:numMatSamplesExternal_int
    for i=1:length( perturbAmount_vec )
        yvec_lasso( i ) = min( min( numMistakes_lasso(i,:,:,ms) ) );
        yvec_inv( i ) = min( min( numMistakes_inv(i,:,:,ms) ) );
    end
    
    fig = plot( perturbAmount_vec, yvec_lasso );
    legend( 'mistakes across perturbations' );
    xlabel( 'perturbation' ); ylabel( 'mistakes' );
    title('min(lamdasxrepititions), avgMat, lasso recovery smethod');
    saveas( fig, [ fname_pert '/lasso_sample_' num2str(ms) ] );

    fig = plot( perturbAmount_vec, yvec_inv );
    legend( 'mistakes across perturbations' );
    title('min(thresholdxrepititions), avgMat, linear recovery method');
    xlabel( 'perturbation' ); ylabel( 'mistakes' );
    saveas( fig, [ fname_pert '/inv_sample_' num2str(ms) ] );
end
