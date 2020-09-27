setup_mapper;
setup_encoder;

%% Simulation parameters.
sim_N = 400;
N = 4096;
SNR_arr = [0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10, 12.5, 15]; % target SNR.
Ps = 1;

%% Start simulation.
SNRs_abs = 10.^(SNR_arr/10);
sigma_arr = sqrt(Ps./SNRs_abs);
err_bit_cnt_after_coding = zeros(length(sigma_arr), 1);
err_bit_cnt_before_coding = zeros(length(sigma_arr), 1);
soft_decode = true;

for sigma_iter = 1:length(sigma_arr)
    sigma = sigma_arr(sigma_iter);
    
    for sim_iter = 1:sim_N
        random_bits = (rand([1, N])>0.5);
        encoded_bits = conv_encode(random_bits, conv_encoder_conf);

        syms = bit_mapping(encoded_bits, mapping_conf);
        ch = ch_realization(length(syms), ch_conf);
        syms_with_noise = syms .* ch + (get_cgaussian(sigma, length(syms))).';
        pred_bits = bit_demapping(syms_with_noise, length(encoded_bits), mapping_conf, ch);
        
        if ~soft_decode
            err_bit_cnt_before_coding(sigma_iter) = err_bit_cnt_before_coding(sigma_iter) + ...
                sum(xor(pred_bits, encoded_bits));
        end
        
        %% Decode.
        decoded_bits = conv_decode(pred_bits, conv_encoder_conf, soft_decode);

        %% Find all the errors.
        err_bit_cnt_after_coding(sigma_iter) = err_bit_cnt_after_coding(sigma_iter) +...
            sum(xor(random_bits, decoded_bits));
    end
    disp([num2str(sigma_iter),'/', num2str(length(sigma_arr)),' Complete for SNR=',...
        num2str(SNR_arr(sigma_iter)), 'dB']);
    if ~soft_decode
        disp(['Log BER before encoding: ', ...
            num2str(log10(sum(err_bit_cnt_before_coding(sigma_iter))/(length(encoded_bits)*sim_N)))]);
    end
    disp(['Log BER after encoding: ', ...
        num2str(log10((sum(err_bit_cnt_after_coding(sigma_iter))/(N*sim_N))))]);
end


%% Display!!
figure(1);
hold on;
set(gca, 'yscale', 'log');
plot(SNR_arr, ((err_bit_cnt_after_coding)/(N*sim_N)).');
plot(SNR_arr, (err_bit_cnt_before_coding/(length(encoded_bits)*sim_N)).');
legend('BER_ conv', 'BER_ ch');
title('BER-SNR Curve');
xlabel('SNR(dB)');
ylabel('BER');
grid on;

disp(['b=', num2str(ch_conf.b), ', rho=',num2str(ch_conf.rho)]);

