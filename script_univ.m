% options
linear = false;
use_lasso_Nmat = false;
n = 10;
A_limiter = 0.9;
numMatrix_samples = 30;
maxIterations = 1000;
itDiff_threshold = 1E-4;
perturb_amount = 0.1;
num_samples = 1;
mistake_threshold = 1E-1;
sigmas = 0:0.025:0.1;
lambdas = 10.^(-6:1:2);
if linear
    func_type = 'linear';
else
    func_type = 'nonlinear';
end

% algorithm starts here
[ numMistakes_avg_Msamples, numMistakes_devs_Msamples ] = sample_numMistakes_nmatrices( n, ...
    itDiff_threshold, mistake_threshold, maxIterations, numMatrix_samples, ...
    perturb_amount, num_samples, lambdas, sigmas, linear, use_lasso_Nmat, ...
    A_limiter )

% output

plot_x = (ones( 1 + length(lambdas), 1 ) * sigmas)';
plot_y = numMistakes_avg_Msamples;
plot_devs = numMistakes_devs_Msamples;
legends = cell( 1, 1 + length(lambdas) );

legends{1} = 'inv';
for lambda_it = 1:length(lambdas)
    legends{ 1 + lambda_it } = ['lasso: lambda=' num2str( lambdas(lambda_it) ) ];
end

errorbar( plot_x, plot_y, plot_devs, ':o' );
legend( 'avg mistakes' );
xlabel('sigma');
ylabel('recovery mistakes');
legend( legends );
title( {[num2str(n) 'X' num2str(n) ' matrix of a ' func_type ' regulatory network, with '...
    num2str(num_samples) ' samples per perturbation'],...
    ['perturb amount=' num2str(perturb_amount) ', '...
    'A limiter=' num2str(A_limiter) ', and ' num2str(numMatrix_samples) ' matrices sampled']} );
%axis( [ sigmas(1) sigmas(end) 0 10] );
axis('auto');
