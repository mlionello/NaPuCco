function [X, stats ] = create_behav(n_pred, timepts, x_mean, x_std)
% CREATE_BEHAV Create simulated data for N tasks observed T times.
% [X, STATS ] = CREATE_BEHAV(N_PRED, TIMEPTS, X_MEAN, X_STD) returns a
% matrix with TIMEPTS timepoints rows and N_PRED + 1 features columns
% with each features sampled form a normal distribution with mean 
% ranging between X_MEAN(1) and X_MEAN(2) and standard deviation ranging
% between X_STD(1) and X_STD(2).
    X = generate_rnd_samples( x_mean, x_std, n_pred, timepts); % generates a matrix
    % of timepts obesrvations from n_pred normal distribution according to the given intervals 
    
    X = cat( 2, ones( timepts, 1), X); %insert intercept as first of the features
    stats = nan;
end