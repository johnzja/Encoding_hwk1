setup_mapper;

%% Start transfer.
N = 160000;
sigma = 0.4;

random_bits = (rand([1, N])>0.5);
syms = bit_mapping(random_bits, mapping_conf);
ch = ch_realization(length(syms), ch_conf);
syms_with_noise = syms .* ch + (get_cgaussian(sigma, length(syms))).';
pred_bits = bit_demapping(syms_with_noise, N, mapping_conf, ch);

error_pattern = xor(pred_bits, random_bits);
pattern_draw(error_pattern);

err_cnt = sum(error_pattern);
disp(['EBR = ', num2str(err_cnt/N*100), '%']);

function pattern_draw(error_pattern)
% draw the error pattern.
    figure(1);
    n_len = floor(sqrt(length(error_pattern)));
    mat = reshape(error_pattern(1:n_len*n_len), [n_len,n_len]);
    imshow(mat);
end