function ch = ch_realization(len, ch_conf)
% ch_conf (struct):
%   b: uncertainty factor, range:(0, 1).
%   rho: channel memory coefficient, range: (0, 1).
%
    % Step1: Check parameter validity.
    assert(ch_conf.b>=0 && ch_conf.b<=1);
    assert(ch_conf.rho>=0 && ch_conf.rho<=1);
    b = ch_conf.b;
    rho = ch_conf.rho;
    % Step2: Calculate all the z's.
    z = get_cgaussian(1, len);
    
    % Step3: Calculate betas.
    beta = filter(sqrt(1-rho^2), [1, -rho], z);
    
    % Step4: Calculate channels.
    ch = sqrt(1-b^2) + b*beta;
    ch = ch.';
end

