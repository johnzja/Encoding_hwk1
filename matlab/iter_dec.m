function [output] = iter_dec(LR, K, S)
    L = LR(1:32);
    R = LR(33:64);
    output = [mod(R+f(L,K,S),2) L];
end

