function m = RSA_decrypt(encrypted_bits, n, d, original_info_len)
    N_bits_per_in_block = floor(log2(n));
    N_bits_per_out_block = ceil(log2(n));
    N_len = length(encrypted_bits);
    assert(mod(N_len, N_bits_per_out_block) == 0, 'Length of encrypted_bits should be multiples of N_bits_per_block!');
    
    d_bin_form = to_binary(d, ceil(log2(d)));
    
    encrypted_bits = reshape(encrypted_bits, N_bits_per_out_block, []);  % arrange the blocks by column.
    [~, N_blocks] = size(encrypted_bits);
    m = zeros(N_bits_per_in_block, N_blocks);
    for block_iter = 1:N_blocks
        one_block = from_binary(encrypted_bits(:,block_iter).');   % within range [0,2^N_bits_per_block-1].
        m(:,block_iter) = to_binary(fast_power_mod(n, one_block, d_bin_form), N_bits_per_in_block).';
    end
    
    m = logical(reshape(m, 1, []));
    m = m(1:original_info_len);
end