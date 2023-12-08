function x = generate_rnd_samples(x_mu_range, x_std_range, n_col, n_rows)
%GENERATE_RND_SAMPLES this function samples n_rows observations from a number n_col of normal distibution
% with mean and std raning in the given intervals.
    x_mu = generate_rnd_interval(x_mu_range, n_col);
    x_std = generate_rnd_interval(x_std_range, n_col);
    x = x_mu + x_std .* randn(n_rows, n_col);
end
