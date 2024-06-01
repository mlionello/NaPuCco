function combined_pvalue = combining_funct( logp_values, opts)
    arguments
        logp_values
        opts.comb_funct = "fisher"
    end
    combined_pvalue = nan;
    switch opts.comb_funct
        case "fisher"
            combined_sum = logp_values{1};
            for k = 2:length(logp_values)
                combined_sum = combined_sum + logp_values{k};
            end
            combined_pvalue = -2 * combined_sum;

        case "tippet"
            combined_pvalue = min(combined_pvalue);
        case "MudholkarGeorge"
            combined_pvalue = sum(log(1-p_values) - ...
                    log(p_values));
            combined_pvalue = combined_pvalue * ...
                1/pi*sqrt( 3*(5*size(p_values, 1)+4)/...
                size(p_values, 1)/(5*size(p_values, 1)+2) );
        otherwise
            error('combining function not recognized among the ones available')
    end
end