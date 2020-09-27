setup_encoder;

%% Simulation parameters.
sim_N = 2;
N = 4000;

%% Start simulation.
for sim_iter = 1:sim_N
    random_bits = (rand([1, N])>0.5);
    encoded_bits = conv_encode(random_bits, conv_encoder_conf);
    
    %% Add random bit error.
    random_bitflip = (rand([1,length(encoded_bits)])>0.99);
    encoded_bits = xor(encoded_bits, random_bitflip);
    bfr = sum(random_bitflip)/length(encoded_bits);
    disp(['BER in channel: ', num2str(bfr*100), '%'])

    %% Decode.
    decoded_bits = conv_decode(encoded_bits, conv_encoder_conf);

    %% Find all the errors.
    err_bit_cnt = sum(xor(random_bits, decoded_bits));
    disp(['BER after decoding: ', num2str(sum(err_bit_cnt)/N*100),'%']);

end





