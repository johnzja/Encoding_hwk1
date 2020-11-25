function [decrypted_bits] = decrypt(encrypted_bits, key, encrypt_method, original_info_bits)
    if strcmp(encrypt_method, 'DES')
        load('S.mat');
        num_blocks = numel(encrypted_bits)/64;
        encrypted_bits = reshape(encrypted_bits, [64, num_blocks]);
        decrypted_bits = zeros(64,num_blocks);
        for k=1:num_blocks
            decrypted_bits(:,k) = DES_decrypt(encrypted_bits(:,k)',key,S)';
        end
        decrypted_bits = reshape(decrypted_bits, [1,numel(decrypted_bits)]);
        decrypted_bits = decrypted_bits(1:original_info_bits);
    end
    if strcmp(encrypt_method, 'RSA')
        n = key(1);
        d = key(3);
        decrypted_bits = RSA_decrypt(encrypted_bits, n, d, original_info_bits);
    end
end

