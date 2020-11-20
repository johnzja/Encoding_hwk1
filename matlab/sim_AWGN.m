%% Joint-simulation in AWGN waveform channel.
clear;
setup_mapper;
setup_wave;
setup_encoder;
addpath('utils');

%% Simulation parameters.
N_sim = 3000;
N_info_bits = 8000;                                                     % 1kB file to transmit.
Ebn0_arr = [20, 15, 12.5, 10, 8, 6, 5, 4, 3, 2.5, 2, 1.5, 1, 0.5, 0, -0.5, -1, -1.5];     % Eb/n0 in dB.

Es = 0.5;                                   % 1 sym in waveform channel: 0.5 energy.
Eb = Es/mapping_conf.bps * (conv_encoder_conf.n/conv_encoder_conf.k);   % Notice: if hard-decision, Eb is different.
n0_arr = Eb ./ (10.^(Ebn0_arr/10));         % n0 in linear scale.
N_n0s = length(n0_arr);

L_encoded = zeros(N_n0s, 1);
L_transmit_syms = zeros(N_n0s, 1);
err_bit_cnt_soft_decode = zeros(N_n0s, 1);
err_bit_cnt_hard_demap = zeros(N_n0s, 1);
err_bit_cnt_hard_decode = zeros(N_n0s, 1);

disp_waveform = false;
disp_psd = false;

tic;
parfor n0_iter = 1:N_n0s
    % Temporary variables for PSD estimators.
    RF_signals_trans = zeros(N_sim, 64088);
    RF_signals_recv = zeros(N_sim, 64088);
    my_mapping_conf = mapping_conf;
    for sim_iter = 1:N_sim
        random_bits = (rand([1, N_info_bits])>0.5);
        encoded_bits = conv_encode(random_bits, conv_encoder_conf);
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
        err_bit_cnt_soft_decode(n0_iter) = err_bit_cnt_soft_decode(n0_iter) + sum(xor(random_bits, decoded_bits_soft));
        
        %% Hard-demap & Hard-decode
        my_mapping_conf.out = 'hard';
        pred_bits = bit_demapping(syms_recv, L_encoded(n0_iter), my_mapping_conf, ch, ch_conf, 1);
        decoded_bits_hard = fast_conv_decode(pred_bits, conv_encoder_conf, false);
        err_bit_cnt_hard_decode(n0_iter) = err_bit_cnt_hard_decode(n0_iter) + sum(xor(random_bits, decoded_bits_hard));
        
        %% Hard-demap only.
        syms_transmit = bit_mapping(random_bits, my_mapping_conf);
        wave_transmit = syms2waveform(syms_transmit, waveform_conf, disp_waveform);
        [syms_recv, wave_recv] = waveform2syms(wave_transmit, n0_arr(n0_iter), length(syms_transmit), ...
            waveform_conf, disp_waveform);
        my_mapping_conf.out = 'hard';
        ch = ones(1, length(syms_recv));
        pred_bits = bit_demapping(syms_recv, N_info_bits, my_mapping_conf, ch, ch_conf, 1);
        err_bit_cnt_hard_demap(n0_iter) = err_bit_cnt_hard_demap(n0_iter) + sum(xor(random_bits, pred_bits));
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
plot(Ebn0_arr, (err_bit_cnt_soft_decode/(N_info_bits*N_sim)).','r-*');
plot(Ebn0_arr, (err_bit_cnt_hard_decode/(N_info_bits*N_sim)).','b-*');
plot(Ebn0_arr-10*log10(2), (err_bit_cnt_hard_demap/(N_info_bits*N_sim)).','g-*');   % Hard-demapping: Eb/n0-3dB
legend('BER_ conv_ soft_ decode', 'BER_ conv_ hard_ decode', 'BER_ no_ coding');
set(gca, 'yscale', 'log');
title('BER-Eb/n_0 Curve');
xlabel('Eb/n0(dB)');
ylabel('BER');
grid on;

disp(['Time elapsed: ', num2str(time_elapsed), 's for ', num2str(N_n0s*N_sim), ...
    ' waveform channel simulations']);

%% SAVE
save(['data/sim_AWGN_soft_', strrep(datestr(datetime), ':', '_'), '.mat']);
