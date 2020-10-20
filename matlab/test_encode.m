clear all;
setup_encoder;
setup_mapper;
conv_encoder_trailing=false;
%% Simulation parameters.
sim_N = 10;
N = 40960;

sigma = 0.3;
mapping_conf.interleave=false;
ch_conf.b = 0;
ch_conf.rho = 0.95;

%% Start simulation.
for sim_iter = 1:sim_N
    random_bits = (rand([1, N])>0.5);
    encoded_bits = conv_encode(random_bits, conv_encoder_conf);
    
    syms = bit_mapping(encoded_bits, mapping_conf);
    ch = ch_realization(length(syms), ch_conf);
    syms_with_noise = syms .* ch + (get_cgaussian(sigma, length(syms))).';
    
    mapping_conf.out = 'L2'; 
    pred_probs = bit_demapping(syms_with_noise, length(encoded_bits), mapping_conf, ch, ch_conf, sigma);
    mapping_conf.out = 'hard'; 
    pred_bits = bit_demapping(syms_with_noise, length(encoded_bits), mapping_conf, ch, ch_conf, sigma);
    %% Decode.
    hard_decoded_bits = fast_conv_decode(pred_bits, conv_encoder_conf, false);
    soft_decoded_bits = fast_conv_decode(pred_probs, conv_encoder_conf, true);

    %% Find all the errors.
    err_bit_soft_cnt = sum(xor(random_bits, soft_decoded_bits));
    err_bit_hard_cnt = sum(xor(random_bits, hard_decoded_bits));
    disp(['BER after decoding hard: ', num2str(sum(err_bit_hard_cnt)/N*100),'%']);
    disp(['BER after decoding soft: ', num2str(sum(err_bit_soft_cnt)/N*100),'%']);

end





