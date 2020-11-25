function [LR] = iter_dec(LR, K, S)
    L = LR(1:32);
    R = LR(33:64);
    LR(33:64) = L;
    LR(1:32) = mod(R+f(L,K,S),2);
end

