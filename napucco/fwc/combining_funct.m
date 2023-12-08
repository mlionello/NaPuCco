function combined_pvalue = combining_funct( p_values, opts)
    %NEED TO WRITE A TEST FUNCTION
    arguments
        p_values (:,:,:) double
        opts.comb_funct string = "fisher"
    end
    combined_pvalue = nan;
    %list_of_methods = ["fisher", "tippet", "MudholkarGeorge"];
    switch opts.comb_funct
        case "fisher"
            combined_pvalue = - 2*sum(log(p_values));
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