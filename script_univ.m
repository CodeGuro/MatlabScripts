% options
linear = false; % flag to determine whether script uses linear or nonlinear network
use_lasso_Nmat = true; % flag to determine whether the quadratic method is used at lasso
n = 5; % size of the matrix (n by n)
A_limiter = 0.8; % limits the generation of nonzeros for matA ([0,1] inclusive)
numMatrix_samples = 20; % number of unique matrices sampled (new matA generated each time)
maxIterations = 1000; % used for nonlinear networks for finding steady states, warning is given in output if this is exceeded
itDiff_threshold = 1E-4; % used for nonlinear networks for finding steady states. (threshold for differences between previous iteration and current one)
perturb_amount = 1; % perturbation amount (applies to both linear & nonlinear)
num_samples = 20; % number of samples taken per sigma (values of noise may differ when sigma > 0) (linear & nonlinear)
mistake_threshold = 1E-1; % values below this in the matrix recovered from matdeltaX are assumed to be 0 (linear & nonlinear)
sigmas = 0:0.03:0.5; % sigmas associated with the level of noise. Plotted in x-dimension
lambdas = [2E-4 0.02 0.06 0.1 1 2E10 ]; % associated with lasso. Change these freely. Use empty vector if you don't want lasso plotted
if linear
    func_type = 'linear';
else
    func_type = 'nonlinear';
end

% algorithm starts here
[ numMisakes_avg_Msamples, numMisakes_devs_Msamples ] = sample_numMistakes_nmatrices( n, ...
    itDiff_threshold, mistake_threshold, maxIterations, numMatrix_samples, ...
    perturb_amount, num_samples, lambdas, sigmas, linear, use_lasso_Nmat, ...
    A_limiter );

% output

plot_x = (ones( 1 + length(lambdas), 1 ) * sigmas)';
plot_y = numMisakes_avg_Msamples;
plot_devs = numMisakes_devs_Msamples;
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
