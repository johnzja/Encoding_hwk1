function c = RSA_encrypt(info_bits, n, e)
% input info_bits: row vector of bits.
% input n = p*q
% input e: gcd(e,phi)=1.

% Typical (n,e,d) for test = (988027, 283, 449467).

    %% Step1: Group the input bits.
    N_bits_per_in_block = floor(log2(n));
    N_bits_per_out_block = ceil(log2(n));
    % zero-padding if necessary.
    num_comp_bits = mod(numel(info_bits),N_bits_per_in_block);
    if num_comp_bits
        info_bits = [info_bits zeros(1,N_bits_per_in_block - num_comp_bits)];
    end
    
    %% Step2: Encrypt.
    e_bin_form = to_binary(e, ceil(log2(e)));
    info_bits = reshape(info_bits, N_bits_per_in_block, []);  % arrange the blocks by column.
    [~, N_blocks] = size(info_bits);
    c = zeros(N_bits_per_out_block, N_blocks);
    for block_iter = 1:N_blocks
        one_block = from_binary(info_bits(:,block_iter).');   % within range [0,2^N_bits_per_block-1].
        c(:,block_iter) = to_binary(fast_power_mod(n, one_block, e_bin_form), N_bits_per_out_block).';
    end
    
    %% Convert to row vectors.
    c = logical(reshape(c, 1, []));
end