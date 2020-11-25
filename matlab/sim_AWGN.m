%% Joint-simulation in AWGN waveform channel.
clear;
setup_mapper;
setup_wave;
setup_encoder;
setup_encrypter;
addpath('utils');

%% Simulation parameters.
N_sim = 3000;
N_info_bits = 8000;                                                     % 1kB file to transmit.
Ebn0_arr = [20, 15, 12.5, 10, 8, 6, 5, 4, 3, 2.5, 2, 1.5, 1, 0.5, 0, -0.5, -1, -1.5];     % Eb/n0 in dB.

Es = 0.5;                                   % 1 sym in waveform channel: 0.5 energy.
Eb = Es/mapping_conf.bps * (conv_encoder_conf.n/conv_encoder_conf.k);   % Notice: if hard-decision, Eb is different.
n0_arr = Eb ./ (10.^(Ebn0_arr/10));         % n0 in linear scale.
if encrypter.enable == true && strcmp(encrypter.method, 'RSA')
    n0_arr = n0_arr * 8000 / 8440;
end
N_n0s = length(n0_arr);

L_encoded = zeros(N_n0s, 1);
L_transmit_syms = zeros(N_n0s, 1);
err_bit_cnt_soft_decode = zeros(N_n0s, 1);
err_bit_cnt_hard_demap = zeros(N_n0s, 1);
err_bit_cnt_hard_decode = zeros(N_n0s, 1);
err_bit_cnt_soft_decode_decrypted = zeros(N_n0s, 1);
err_bit_cnt_hard_decode_decrypted = zeros(N_n0s, 1);
err_block_cnt_soft_decode_decrypted = zeros(N_n0s, 1);
err_block_cnt_hard_decode_decrypted = zeros(N_n0s, 1);
err_block_cnt_soft_decode = zeros(N_n0s, 1);
err_block_cnt_hard_decode = zeros(N_n0s, 1);

disp_waveform = false;
disp_psd = false;

tic;
parfor n0_iter = 1:N_n0s
    % Temporary variables for PSD estimators.
    if encrypter.enable == true && strcmp(encrypter.method, 'RSA')
        RF_signals_trans = zeros(N_sim, 67608);
        RF_signals_recv = zeros(N_sim, 67608);
    else
        RF_signals_trans = zeros(N_sim, 64088);
        RF_signals_recv = zeros(N_sim, 64088);        
    end
    my_mapping_conf = mapping_conf;
    for sim_iter = 1:N_sim
        random_bits = (rand([1, N_info_bits])>0.5); % Generate info bits. 
        if encrypter.enable == true
            bits_before_encode = encrypt(random_bits,encrypter.method,encrypter.key); % Encryption. 
        else
            bits_before_encode = random_bits;
        end
        encoded_bits = conv_encode(bits_before_encode, conv_encoder_conf);
        L_encoded(n0_iter) = length(encoded_bits);
        syms_transmit = bit_mapping(encoded_bits, my_mapping_conf);  
        L_transmit_syms(n0_iter) = length(syms_transmit);
        wave_transmit = syms2waveform(syms_transmit, waveform_conf, disp_waveform);
        RF_signals_trans(sim_iter, :) = wave_transmit;  % RECORD
        
        [syms_recv, wave_recv] = waveform2syms(wave_transmit, n0_arr(n0_iter), length(syms_transmit), ...
            waveform_conf, disp_waveform);
        RF_signals_recv(sim_iter, :) = wave_recv;       % RECORD
        
        %% Soft-demap & Soft-decode
        ch = ones(1, length(syms_recv));
        my_mapping_conf.out = 'L2';
        pred_L2 = bit_demapping(syms_recv, L_encoded(n0_iter), my_mapping_conf, ch, ch_conf, 1);
        decoded_bits_soft = fast_conv_decode(pred_L2, conv_encoder_conf, true); % Perform soft-decode.
        err_bit_cnt_soft_decode(n0_iter) = err_bit_cnt_soft_decode(n0_iter) + sum(xor(bits_before_encode, decoded_bits_soft));
        
        %% Hard-demap & Hard-decode
        my_mapping_conf.out = 'hard';
        pred_bits = bit_demapping(syms_recv, L_encoded(n0_iter), my_mapping_conf, ch, ch_conf, 1);
        decoded_bits_hard = fast_conv_decode(pred_bits, conv_encoder_conf, false);
        err_bit_cnt_hard_decode(n0_iter) = err_bit_cnt_hard_decode(n0_iter) + sum(xor(bits_before_encode, decoded_bits_hard));
        
        %% Hard-demap only.
        syms_transmit = bit_mapping(random_bits, my_mapping_conf);
        wave_transmit = syms2waveform(syms_transmit, waveform_conf, disp_waveform);
        [syms_recv, wave_recv] = waveform2syms(wave_transmit, n0_arr(n0_iter), length(syms_transmit), ...
            waveform_conf, disp_waveform);
        my_mapping_conf.out = 'hard';
        ch = ones(1, length(syms_recv));
        pred_bits = bit_demapping(syms_recv, N_info_bits, my_mapping_conf, ch, ch_conf, 1);
        err_bit_cnt_hard_demap(n0_iter) = err_bit_cnt_hard_demap(n0_iter) + sum(xor(random_bits, pred_bits));
        blk_size = encrypter.blk_size_encrypted_bits;
        for i = 1:blk_size:length(decoded_bits_soft)
            blk_start = i;
            blk_end = min(i+blk_size-1, length(decoded_bits_soft));
            err_block_cnt_soft_decode(n0_iter) = err_block_cnt_soft_decode(n0_iter) +...
                ~isequal(bits_before_encode(blk_start:blk_end), decoded_bits_soft(blk_start:blk_end));  
            err_block_cnt_hard_decode(n0_iter) = err_block_cnt_hard_decode(n0_iter) +...
                ~isequal(bits_before_encode(blk_start:blk_end), decoded_bits_hard(blk_start:blk_end));  
        end
        
        %% Decrypt
        if encrypter.enable == true
            decoded_bits_soft_decrypted = decrypt(decoded_bits_soft,encrypter.key,encrypter.method,N_info_bits);
            decoded_bits_hard_decrypted = decrypt(decoded_bits_hard,encrypter.key,encrypter.method,N_info_bits);
            err_bit_cnt_soft_decode_decrypted(n0_iter) = ...
                err_bit_cnt_soft_decode_decrypted(n0_iter) + sum(xor(random_bits, decoded_bits_soft_decrypted));
            err_bit_cnt_hard_decode_decrypted(n0_iter) = ...
                err_bit_cnt_hard_decode_decrypted(n0_iter) + sum(xor(random_bits, decoded_bits_hard_decrypted));
            blk_size = encrypter.blk_size_info_bits;
            for i = 1:blk_size:length(random_bits)
                blk_start = i;
                blk_end = min(i+blk_size-1, length(random_bits));
                err_block_cnt_soft_decode_decrypted(n0_iter) = err_block_cnt_soft_decode_decrypted(n0_iter) +...
                    ~isequal(random_bits(blk_start:blk_end), decoded_bits_soft_decrypted(blk_start:blk_end));  
                err_block_cnt_hard_decode_decrypted(n0_iter) = err_block_cnt_hard_decode_decrypted(n0_iter) +...
                    ~isequal(random_bits(blk_start:blk_end), decoded_bits_hard_decrypted(blk_start:blk_end)); 
            end
        end
        
       % Display running info.
        if mod(sim_iter, floor(N_sim/10))==0
            disp(['SOFT Eb/n0=', num2str(Ebn0_arr(n0_iter)),'dB: ',num2str(sim_iter/N_sim*100),'% complete']);
        end
    end
    disp([num2str(n0_iter),'/', num2str(N_n0s),' Complete for Eb/n0=',...
        num2str(Ebn0_arr(n0_iter)), 'dB']);
    disp(['Log BER after encoding (soft demapping): ', ...
        num2str(log10(err_bit_cnt_soft_decode(n0_iter)/(N_info_bits*N_sim)))]);
    disp(['Log BER after encoding (hard demapping): ', ...
        num2str(log10(err_bit_cnt_hard_decode(n0_iter)/(N_info_bits*N_sim)))]);
    disp(['Log BER without encoding (hard demapping): ', ...
        num2str(log10(err_bit_cnt_hard_demap(n0_iter)/(N_info_bits*N_sim)))]);
    disp(['Log BER with encoding and encryption (soft demapping): ', ...
        num2str(log10(err_bit_cnt_soft_decode_decrypted(n0_iter)/(N_info_bits*N_sim)))]);

    %% Draw PSD.
    if disp_psd
        figure;
        psd_est(RF_signals_trans, waveform_conf.oversample_rate * L_transmit_syms(n0_iter), waveform_conf.fs, true);
        hold on;
        psd_est(RF_signals_recv, waveform_conf.oversample_rate * L_transmit_syms(n0_iter), waveform_conf.fs, true);
        legend('trans PSD', 'recv PSD');
    end
end
time_elapsed = toc;

%% Disp
figure;
hold on;
plot(Ebn0_arr, (err_bit_cnt_soft_decode/(N_info_bits*N_sim)).', '-*');
plot(Ebn0_arr, (err_bit_cnt_hard_decode/(N_info_bits*N_sim)).', '-*');
Ebn0_arr_hard_demapping = Ebn0_arr - 10*log10(2);
Ebn0_arr_hard_demapping = Ebn0_arr_hard_demapping(err_bit_cnt_hard_demap~=0);
err_bit_cnt_hard_demap=err_bit_cnt_hard_demap(err_bit_cnt_hard_demap~=0);
plot(Ebn0_arr_hard_demapping, (err_bit_cnt_hard_demap/(N_info_bits*N_sim)).','-*');   % Hard-demapping: Eb/n0-3dB
plot(Ebn0_arr, (err_bit_cnt_soft_decode_decrypted/(N_info_bits*N_sim)).', '-*');
BER_theoretical = qfunc(sqrt(2*10.^(Ebn0_arr_hard_demapping/10)));
plot(Ebn0_arr_hard_demapping, BER_theoretical, '-*');
legend('BER_ conv_ soft_ decode', 'BER_ conv_ hard_ decode', 'BER_ no_ coding',  'BER_ soft_ decode_ with_ encryption', 'BER theoretical');
set(gca, 'yscale', 'log');
title('BER-Eb/n_0 Curve');
xlabel('Eb/n0(dB)');
ylabel('BER');
grid on;

figure;
hold on;
plot(Ebn0_arr, (err_block_cnt_soft_decode_decrypted/((ceil(N_info_bits)/64)*N_sim)).', '-+');
plot(Ebn0_arr, (err_block_cnt_hard_decode_decrypted/((ceil(N_info_bits)/64)*N_sim)).', '-*');
plot(Ebn0_arr, (err_block_cnt_soft_decode/((ceil(N_info_bits)/64)*N_sim)).', '-x');
plot(Ebn0_arr, (err_block_cnt_hard_decode/((ceil(N_info_bits)/64)*N_sim)).','-o');
legend('BLER_ conv_ soft_ decode_ RSA', 'BLER_ conv_ hard_ decode_ RSA', 'BLER_ conv_ soft_ decode',  'BLER_ conv_ hard_ decode');
set(gca, 'yscale', 'log');
title('BLER-Eb/n_0 Curve');
xlabel('Eb/n0(dB)');
ylabel('BLER');
grid on;

disp(['Time elapsed: ', num2str(time_elapsed), 's for ', num2str(N_n0s*N_sim), ...
    ' waveform channel simulations']);

%% SAVE files!
save(['data/sim_AWGN_soft_', strrep(datestr(datetime), ':', '_'), '.mat']);
