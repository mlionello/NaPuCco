function [y, B, R2_current, stats] = createSingVoxBehav( ...
    X, R2_target, R2_tollerance, vx_std_range, logger, opts)
    arguments
    X double;
    R2_target double;
    R2_tollerance double = 0.01;
    vx_std_range double = 1;
    logger = nan;
    opts.verbose int32 = 0;
    opts.n_iter int32 = 100;
    opts.initial_beta_range double = [-0.07 0.005];
    opts.step_in_beta double = 0.005;
    end
    stats = nan;

    [n_timepoints, n_predictors] = size(X);
    beta_range = opts.initial_beta_range;
    
    R2_current = 0;
    B = nan;
    k = 0;
    
    while abs(R2_target - R2_current) > R2_tollerance || k == 0
        k = k + 1;

        betas = generate_rnd_interval(beta_range, n_predictors).';
        for i = 1: opts.n_iter
            
            y = (vx_std_range .* randn(n_timepoints, 1)) + X*betas;

%             [B, ~, ~, ~, STATS] = regress(y, X);
%             R2_current = STATS(1);

            B = X\y;
            yhat = X*B;
            ssres = sum((y-yhat).^2);
            sstot = sum((y-mean(y)).^2);

            R2_current = 1-(ssres./sstot);
            
            if abs(R2_target - R2_current) < R2_tollerance
                break
            end
        end
                
        if opts.verbose == 1 
            fprintf('R^2==%.3f distance from required = %.3f\n', R2_current, abs(R2_target - R2_current))
        end
        if mod(k, 100)==0
            R2_tollerance = R2_tollerance*1.7;
        end        
        if R2_current >= min((R2_target + R2_tollerance)*1.3, 0.99)
            beta_range = opts.initial_beta_range;
            continue
        end
        beta_range = beta_range + opts.step_in_beta;
    end
    if opts.verbose > 0
        fprintf('R^2==%.3f distance from required = %.3f\n', R2_current, abs(R2_target - R2_current))
    end
    stats = [R2_current, abs(R2_target - R2_current), R2_tollerance,...
        beta_range, k, i, mean(B), std(B), kurtosis(B), max(B), min(B)];
    logger.tolog(string(stats));
    onlydatalog = fopen(fullfile( ...
        dir(logger.logfile).folder, 'datalog.csv'), 'a');
    fprintf(onlydatalog, '%.5f, ', R2_target, R2_current, R2_tollerance, ...
        beta_range, k, i, mean(B), std(B), kurtosis(B), max(B), min(B), vx_std_range);
    fprintf(onlydatalog, '\n');
    fclose(onlydatalog);
end