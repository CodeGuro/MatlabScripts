load( 'ws_data_2.mat' );

fname_pert = 'perturbations';
fname_reps = 'repititions';

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
    
    fname_rep_thresh = [ fname_reps '/bisect/mat' num2str(ms) ];
    mkdir( fname_rep_thresh );
    
    for perturbidx = 1:length( perturbAmount_vec )
        fname_rep_thresh = [fname_rep_thresh '/perturb_' num2str(perturbidx)];
        mkdir(fname_rep_thresh);
        for threshidx = 1:length(mistakeThresh_real_vec)
            yvec_inv_thresh = numMistakes_inv( perturbidx, :, threshidx, ms );
            fig = plot( numsamples_vec, yvec_inv_thresh );
            title( [ 'threshold = ' num2str( mistakeThresh_real_vec(threshidx) )...
                ' at perturb=' num2str(perturbAmount_vec(perturbidx)) ] );
            xlabel( 'repetitions' );
            ylabel( 'recov mistakes' );
            saveas( fig, [ fname_rep_thresh '/thresh_' num2str( threshidx ) ] );
        end
        
        for lambdaidx = 1:length(lambdas_vecreal)
            yvec_lasso_lambda = numMistakes_lasso( perturbidx, :, lambdaidx, ms );
            fig = plot( numsamples_vec, yvec_lasso_lambda );
            title( [ 'lambda = ' num2str( lambdas_vecreal(lambdaidx) )...
                ' at perturb=' num2str(perturbAmount_vec(perturbidx)) ] );
            xlabel( 'repetitions' );
            ylabel( 'recov mistakes' );
            saveas( fig, [ fname_rep_thresh '/lambda_' num2str( lambdaidx ) ] );
        end
    end
    
    
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
    xlabel( 'samples' ); ylabel( 'mistakes' );
    title('min(lamdasxrepititions), avgMat, lasso recovery smethod');
    saveas( fig, [ fname_pert '/lasso_sample_' num2str(ms) ] );

    fig = plot( perturbAmount_vec, yvec_inv );
    legend( 'mistakes across perturbations' );
    title('min(thresholdxrepititions), avgMat, linear recovery method');
    xlabel( 'samples' ); ylabel( 'mistakes' );
    saveas( fig, [ fname_pert '/inv_sample_' num2str(ms) ] );
end
