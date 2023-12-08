function out = cubicrt(inputArg1)
    out = power(double(inputArg1), 1/3);
    for i=1:length(inputArg1)
        if int32(out(i))^3==inputArg1(i)
            out(i) = int32(out(i));
        end
        out(i) = floor(out(i));
    end
end