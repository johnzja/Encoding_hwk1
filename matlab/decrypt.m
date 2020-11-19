function [decrypted_bits] = decrypt(encrypted_bits, key, encrypt_method)
    if strcmp(encrypt_method, 'DES')
        load('S.mat');
        num_blocks = numel(encrypted_bits)/64;
        encrypted_bits = reshape(encrypted_bits, [64, num_blocks]);
        decrypted_bits = zeros(64,num_blocks);
        for k=1:num_blocks
            decrypted_bits(:,k) = DES_decrypt(encrypted_bits(:,k)',key,S)';
        end
        decrypted_bits = reshape(decrypted_bits, [1,numel(decrypted_bits)]);
    end
end

