function p = fast_power_mod(n, a, alpha_binform)
% n: modulo
% a: base
% alpha_binform: alpha in binary form.
    alpha_bits = length(alpha_binform);
    p = 1;
    apower = a;
    for k = alpha_bits:-1:1
        if alpha_binform(k)
            p = mod(p * apower, n);
        end
        apower = mod(apower * apower, n);
    end
end