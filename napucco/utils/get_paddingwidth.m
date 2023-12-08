function m = get_paddingwidth(n, b)
    m = roots([8, double(12*n), double(6*n*n), -double(b)]);
    m = real(m(end));
end
