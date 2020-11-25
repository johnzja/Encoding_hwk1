function [encrypted_bits] = encrypt(info_bits, encrypt_method, key)
    if strcmp(encrypt_method, 'DES')
        load('S.mat');
        num_comp_bits = mod(numel(info_bits),64);
        if num_comp_bits
            info_bits = [info_bits zeros(1,64 - num_comp_bits)];
        end
        num_blocks = numel(info_bits)/64;
        info_bits = reshape(info_bits, [64, num_blocks]);
        encrypted_bits = zeros(64,num_blocks);
        for k=1:num_blocks
            encrypted_bits(:,k) = DES_encrypt(info_bits(:,k)',key,S)';
        end
        encrypted_bits = reshape(encrypted_bits, [1,numel(encrypted_bits)]);
    end
    if strcmp(encrypt_method, 'RSA')
        n = key(1);
        e = key(2);
        encrypted_bits = RSA_encrypt(info_bits, n, e);
    end
end

