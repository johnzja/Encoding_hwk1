clear all;
load('data/rcf.mat');

disp_flag = false;
%% Generate complex symbols with energy 1.
N_syms = 100;
N_sim = 1000;

oversample_rate = 8;
fs = 16000;     % sample rate
fc = 1850;      % carrier freq = 1850Hz
sigma = 0.1;
len_signal = N_syms * oversample_rate + 2*Group_delay;

transmit_delta_sequence = zeros(len_signal, 1);
transmit_signal_RF= cell(N_sim, 1);

for sim_iter = 1:N_sim
    syms_I = 1-2*(rand([N_syms, 1])>0.5);
    syms_Q = 1-2*(rand([N_syms, 1])>0.5);

    for k = 1:N_syms
        transmit_delta_sequence(1+(k-1)*oversample_rate) = (syms_I(k) + 1j*syms_Q(k))/sqrt(2);
    end
    transmit_signal_baseband = filter(g_arr, [1], transmit_delta_sequence);  % generate I/Q signals  && put into physical AWGN channel.
    if disp_flag
        figure(1);
        subplot(2,1,1);
        n_arr = (0:length(transmit_delta_sequence)-1).';
        stem(n_arr, real(transmit_delta_sequence));
        hold on;
        stem(n_arr, real(transmit_signal_baseband));
        title('(Re) delta sequence && x_B_B_,_I(t)');
        legend('transmit delta', 'transmit baseband');
    end
    % up-convert.
    % channel parameters.

    transmit_signal = real(transmit_signal_baseband .* exp(1j*2*pi*((0:len_signal-1).')*fc/fs));    % upconvert.
    transmit_signal_RF{sim_iter} = transmit_signal; % Record RF signal.
    
    if disp_flag
        subplot(2,1,2);
        stem(n_arr, transmit_signal);
        title('(RF) signal transmitted within band 300~3400Hz');
        figure(2);
    end
    
    noises = sigma * randn(size(transmit_signal));
    recv_signal = transmit_signal + noises;                                                         % Add AWGN noise.
    recv_signal_baseband = (2*transmit_signal).*exp(-1j*2*pi*((0:len_signal-1).')*fc/fs);
    recv_noise_baseband = (2*noises).*exp(-1j*2*pi*((0:len_signal-1).')*fc/fs);

    recv_signal_after_MF = filter(g_arr, [1], recv_signal_baseband);                                 % g_arr filter is an LPF.
    recv_noise_after_MF = filter(g_arr, [1], recv_noise_baseband);
    
    if disp_flag
        subplot(2,1,1);
        stem(real(recv_signal_after_MF));
        hold on;
        stem(real(recv_signal_after_MF + recv_noise_after_MF));
        legend('recv BB', 'recv BB with noise');
    end
    disp(['Energy_syms=', num2str(sum(abs(transmit_delta_sequence).^2))]);
    disp(['Energy_transmitted_signal=', num2str(sum(abs(transmit_signal).^2))]);
    disp(['channel Es/n0(dB) = ', num2str(0.5*db(0.5/(2*sigma^2)))]);

    recv_delay = Group_delay*2;
    recv_sample_signal_after_MF = recv_signal_after_MF(recv_delay+1:oversample_rate:end);       % complex
    recv_sample_noise_after_MF = recv_noise_after_MF(recv_delay+1:oversample_rate:end);         % complex
    disp(['Energy_sampled_noise_after_MF=', num2str(sum(abs(recv_sample_noise_after_MF).^2))]);
    disp(['Energy_sampled_signal_after_MF=', num2str(sum(abs(recv_sample_signal_after_MF).^2))]);
    disp(['statistical SNR(dB) = ', num2str(...
        1/2*db(sum(abs(recv_sample_signal_after_MF).^2)/sum(abs(recv_sample_noise_after_MF).^2)))]);
    % Expected: SNR <= Es/n0 + 3dB
end

%% Analyze RF signal spectrum.
% Spectrum estimation.
addpath('utils');
% construct matrix.
sample_len = length(transmit_signal_RF{1});
RF_observation_matrix = zeros(N_sim, sample_len);
for k = 1:N_sim
    RF_observation_matrix(k,:) = (transmit_signal_RF{k}).';
end

psd_est(RF_observation_matrix, oversample_rate * N_syms, fs, true);

