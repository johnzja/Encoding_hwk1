function [seq] = DES_encrypt(seq, K, S)
    Kn = KS(K);
    seq = IP(seq);
    for i=1:16
        seq = iter_enc(seq, Kn(i,:), S);
    end
    seq = invIP(seq);
end

