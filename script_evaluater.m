load( 'ws_data_updated.mat' );

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
saveas( fig, 'folder/lasso_avg' );

fig = plot( numsamples_vec, yvec_inv );
legend( 'mistakes across repetitions' );
title('min(thresholdxperturb), avgMat, linear recovery method');
xlabel( 'samples' ); ylabel( 'mistakes' );
saveas( fig, 'folder/inv_avg' );

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
    saveas( fig, [ 'folder/lasso_sample_' num2str(ms) ] );

    fig = plot( numsamples_vec, yvec_inv );
    legend( 'mistakes across repetitions' );
    title('min(thresholdxperturb), avgMat, linear recovery method');
    xlabel( 'samples' ); ylabel( 'mistakes' );
    saveas( fig, [ 'folder/inv_sample_' num2str(ms) ] );
    
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
saveas( fig, 'folder2/lasso_avg' );

fig = plot( perturbAmount_vec, yvec_inv );
legend( 'mistakes across perturbations' );
title('min(thresholdxrepititions), avgMat, linear recovery method');
xlabel( 'samples' ); ylabel( 'mistakes' );
saveas( fig, 'folder2/inv_avg' );

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
    saveas( fig, [ 'folder2/lasso_sample_' num2str(ms) ] );

    fig = plot( perturbAmount_vec, yvec_inv );
    legend( 'mistakes across perturbations' );
    title('min(thresholdxrepititions), avgMat, linear recovery method');
    xlabel( 'samples' ); ylabel( 'mistakes' );
    saveas( fig, [ 'folder2/inv_sample_' num2str(ms) ] );
    
end