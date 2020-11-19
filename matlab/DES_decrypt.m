function [seq] = DES_decrypt(seq, K, S)
    Kn = KS(K);
    seq = IP(seq);
    for i=1:16
        seq = iter_dec(seq, Kn(17-i,:), S);
    end
    seq = invIP(seq);
end

