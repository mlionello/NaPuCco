function x = generate_rnd_interval(interval, n_points)
    x = interval(1) + (interval(2) - interval(1)) .* rand(1, n_points);
end