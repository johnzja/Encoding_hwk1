function [LR] = iter_enc(LR, K, S)
    L = LR(1:32);
    R = LR(33:64);
    LR(1:32) = R;
    LR(33:64) = mod(L+f(R,K,S),2);
end

