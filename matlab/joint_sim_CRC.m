clear;

setup_mapper;
setup_encoder;

%% Simulation parameters.
N_sim = 1000;
N_info_bits = 4096;
SNR_arr = [0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 12.5, 15, 17.5];   % target SNR.
Ps = 1;
soft_decode = true;
record_csi = false;

%% Start simulation.
SNRs_abs = 10.^(SNR_arr/10);
sigma_arr = sqrt(Ps./SNRs_abs);
N_sigmas = length(sigma_arr);

err_bit_cnt_after_hard_decoding = zeros(N_sigmas, 1);
err_bit_cnt_before_coding = zeros(N_sigmas, 1);
err_box_cnt_crc_hard = zeros(N_sigmas,1);
L_encoded = zeros(N_sigmas, 1);

if soft_decode
    err_bit_cnt_after_soft_decoding = zeros(N_sigmas, 1);
    err_box_cnt_crc_soft = zeros(N_sigmas,1);
end

% Record all the channel CSI and Constellation points.
CSI = cell(N_sigmas, 1);    % each row of each cell element: The ch.
SYMS_TRANSMIT = cell(N_sigmas, 1);
SYMS_RECEIVE = cell(N_sigmas, 1);

tic;
parfor sigma_iter = 1:N_sigmas
    sigma = sigma_arr(sigma_iter);
    
    for sim_iter = 1:N_sim
        random_bits = (rand([1, N_info_bits])>0.5);
        random_bits_with_crc = CRC(random_bits);
        encoded_bits = conv_encode(random_bits_with_crc, conv_encoder_conf);
        
        %% Setup length of encoded bits.
        L_encoded(sigma_iter) = length(encoded_bits);
        
        %% Simulation with channel.
        syms = bit_mapping(encoded_bits, mapping_conf);
        ch = ch_realization(length(syms), ch_conf);
        syms_with_noise = syms .* ch + (get_cgaussian(sigma, length(syms))).';
        if record_csi
            if isempty(CSI{sigma_iter})
                Ls = length(syms);
                CSI{sigma_iter} = zeros(N_sim, Ls);
                SYMS_TRANSMIT{sigma_iter} = zeros(N_sim, Ls);
                SYMS_RECEIVE{sigma_iter} = zeros(N_sim, Ls);
            end
            CSI{sigma_iter}(sim_iter, :) = ch;
            SYMS_TRANSMIT{sigma_iter}(sim_iter, :) = syms;
            SYMS_RECEIVE{sigma_iter}(sim_iter, :) = syms_with_noise;
        end
        
        my_mapping_conf = mapping_conf;     % Copy mapping_conf as a temp-variable to enable parfor.
        if soft_decode
            my_mapping_conf.out = 'L2';     % Switch to soft-demapping.
            pred_probs = bit_demapping(syms_with_noise, length(encoded_bits), my_mapping_conf, ch, ch_conf, sigma);
        end
        my_mapping_conf.out = 'hard';       % Switch to hard-demapping.
        pred_bits = bit_demapping(syms_with_noise, length(encoded_bits), my_mapping_conf, ch, ch_conf, sigma);
        
        err_bit_cnt_before_coding(sigma_iter) = err_bit_cnt_before_coding(sigma_iter) + ...
            sum(xor(pred_bits, encoded_bits));
        
        %% Decode.
        if soft_decode
            soft_decoded_bits_with_crc = fast_conv_decode(pred_probs, conv_encoder_conf, soft_decode);
            soft_decoded_validation = deCRC(soft_decoded_bits_with_crc);
        end
        hard_decoded_bits_with_crc = fast_conv_decode(pred_bits, conv_encoder_conf, false); % Perform hard-decode.
        hard_decoded_validation = deCRC(hard_decoded_bits_with_crc);
        
        %% Find all the errors.       
        err_bit_cnt_after_hard_decoding(sigma_iter) = err_bit_cnt_after_hard_decoding(sigma_iter) + ...
        sum(xor(random_bits_with_crc, hard_decoded_bits_with_crc));
        if soft_decode
            err_bit_cnt_after_soft_decoding(sigma_iter) = err_bit_cnt_after_soft_decoding(sigma_iter) + ...
                sum(xor(random_bits_with_crc, soft_decoded_bits_with_crc));
            
            for k=1:length(soft_decoded_validation(:,1))
                if(any(soft_decoded_validation(k, :)))
                    err_box_cnt_crc_soft(sigma_iter)=err_box_cnt_crc_soft(sigma_iter)+1;
                end
            end
        end
        for k=1:length(hard_decoded_validation(:, 1))
            if(any(hard_decoded_validation(k, :)))
                err_box_cnt_crc_hard(sigma_iter)=err_box_cnt_crc_hard(sigma_iter)+1;
            end
        end
        
       % Display running info.
        if mod(sim_iter, floor(N_sim/10))==0
            disp(['SNR=', num2str(SNR_arr(sigma_iter)),': ',num2str(sim_iter/N_sim*100),'% complete']);
        end
        
    end
    
    % Simulation loop ends here.
    err_box_cnt_crc_hard(sigma_iter)=err_box_cnt_crc_hard(sigma_iter)/...
        (length(hard_decoded_validation)*N_sim);
    if soft_decode
        err_box_cnt_crc_soft(sigma_iter)=err_box_cnt_crc_soft(sigma_iter)/...
            (length(hard_decoded_validation)*N_sim);
    end
    
    disp([num2str(sigma_iter),'/', num2str(N_sigmas),' Complete for SNR=',...
        num2str(SNR_arr(sigma_iter)), 'dB']);

    disp(['Log BER before encoding (hard demapping): ', ...
        num2str(log10(err_bit_cnt_before_coding(sigma_iter)/(L_encoded(sigma_iter)*N_sim)))]);

    disp(['Log BER after encoding (hard demapping): ', ...
        num2str(log10(err_bit_cnt_after_hard_decoding(sigma_iter)/(N_info_bits*N_sim)))]);
    disp(['Log BLER after encoding (hard demapping): ', ...
        num2str(log10(err_box_cnt_crc_hard(sigma_iter)))]);
    
    if soft_decode
        disp(['Log BER after encoding (soft demapping): ', ...
            num2str(log10(err_bit_cnt_after_soft_decoding(sigma_iter)/(N_info_bits*N_sim)))]);
        disp(['Log BLER after encoding (soft demapping): ', ...
            num2str(log10(err_box_cnt_crc_soft(sigma_iter)))]);
    end
    
end
time_elapsed = toc;
assert(~any(diff(L_encoded)), 'error in length of encoded bits!');


%% Count BLER.
figure(1);
hold on;
set(gca, 'yscale', 'log');
plot(SNR_arr,err_box_cnt_crc_hard.');
if soft_decode
    plot(SNR_arr,err_box_cnt_crc_soft.');
end

if soft_decode
    legend('BLER_ conv_ hard', 'BLER_ conv_ soft');
else
    legend('BLER_ conv_ hard');
end
title('BLER-SNR Curve');

ylabel('BLER');
xlabel('SNR_(_d_B_)');
grid on;

%% Display!!
figure(2);
hold on;
plot(SNR_arr, ((err_bit_cnt_after_hard_decoding)/(N_info_bits*N_sim)).');   % Hard decode is always done.
if soft_decode
    plot(SNR_arr, ((err_bit_cnt_after_soft_decoding)/(N_info_bits*N_sim)).');
end

plot(SNR_arr, (err_bit_cnt_before_coding/(L_encoded(1)*N_sim)).');

if soft_decode
    legend('BER_ conv_ hard', 'BER_ conv_ soft', 'BER_ ch');
else
    legend('BER_ conv_ hard', 'BER_ ch');
end
set(gca, 'yscale', 'log');

title('BER-SNR Curve');
xlabel('SNR(dB)');
ylabel('BER');
grid on;

disp(['Time elapsed: ', num2str(time_elapsed), 's for ', num2str(N_sigmas*N_sim), ' channel simulations']);
disp(['b=', num2str(ch_conf.b), ', rho=',num2str(ch_conf.rho)]);

%% Save variables into files.
save(['data/sim_', strrep(datestr(datetime), ':', '_'), '.mat']);
