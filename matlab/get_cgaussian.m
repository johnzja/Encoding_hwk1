function n=get_cgaussian(sigma, len)
% n=nr+j*ni, E(|n|^2)=sigma^2.
    n = sigma * (randn([len, 1])+1j*randn([len, 1]))/sqrt(2);
end