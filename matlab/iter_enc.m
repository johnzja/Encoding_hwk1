function [output] = iter_enc(LR, K, S)
    L = LR(1:32);
    R = LR(33:64);
    output = [R mod(L+f(R,K,S),2)];
end

