% options
linear = false;
use_lasso_Nmat = false;
nsizes = 3:10;
A_limiter = 0.8;
numMatrix_samples = 12;
maxIterations = 1000;
itDiff_threshold = 1E-4;
perturb_amount = 1;
num_samples = 5;
mistake_threshold = 1E-1;
sigmas = [ 0.03 ];
lambdas = [2E-4 0.02 0.06 0.1 0.2 2E10 ];
if linear
    func_type = 'linear';
else
    func_type = 'nonlinear';
end

% algorithm starts here
numMistakes_avg_nSizes = nan( length(nsizes), 1 + length(lambdas) );
numMistakes_std_nSizes = nan( length(nsizes), 1 + length(lambdas) );
for n_it = 1: length(nsizes)
    n = nsizes(n_it);
    [ numMisakes_avg_Msamples, numMisakes_devs_Msamples ] = sample_numMistakes_nmatrices( n, ...
        itDiff_threshold, mistake_threshold, maxIterations, numMatrix_samples, ...
        perturb_amount, num_samples, lambdas, sigmas, linear, use_lasso_Nmat, ...
        A_limiter );
    numMistakes_avg_nSizes( n_it, : ) = numMisakes_avg_Msamples;
    numMistakes_std_nSizes( n_it, : ) = numMisakes_devs_Msamples;
end

% output
plot_x = (ones( 1 + length(lambdas), 1 ) * nsizes)';
plot_y = numMistakes_avg_nSizes;
plot_devs = numMistakes_std_nSizes;
legends = cell( 1, 1 + length(lambdas) );

legends{1} = 'inv';
for lambda_it = 1:length(lambdas)
    legends{ 1 + lambda_it } = ['lasso: lambda=' num2str( lambdas(lambda_it) ) ];
end

errorbar( plot_x, plot_y, plot_devs, ':o' );
legend( 'avg mistakes' );
xlabel('n size');
ylabel('recovery mistakes');
legend( legends );
title( {[num2str(n) 'X' num2str(n) ' matrix of a ' func_type ' regulatory network, with '...
    num2str(num_samples) ' samples per perturbation'],...
    ['perturb amount=' num2str(perturb_amount) ', sigma=' num2str(sigmas) ', '  ...
    'A limiter=' num2str(A_limiter) ', and ' num2str(numMatrix_samples) ' matrices sampled']} );
%axis( [ sigmas(1) sigmas(end) 0 10] );
axis('auto');
