clear;
load('data/rcf.mat');

% [h, w] = freqz(g_arr, [1], 200);
% figure(1);
% plot(w/pi,abs(h),'.-')
% axis([0 1 -1 2])
% legend('Response');
% ylabel('Magnitude')

%% Generate real symbols with energy 1. (final version: complex signals).
N_syms = 100;
oversample_rate = 4;
len_signal = N_syms * oversample_rate + 40;
transmit_delta_sequence = zeros(len_signal, 1);

% real first.
syms = 1-2*(rand([N_syms, 1])>0.5);
for k = 1:length(syms) 
    transmit_delta_sequence(1+(k-1)*oversample_rate) = syms(k);
end
transmit_signal = filter(g_arr, [1], transmit_delta_sequence);  % Signal put into physical AWGN channel.
figure(1);
subplot(2,1,1);
stem(transmit_delta_sequence);
hold on;
stem(transmit_signal);

sigma = 0.2;
noises = sigma * randn(size(transmit_signal));
recv_signal = transmit_signal + noises;                         % Add AWGN noise.
recv_signal_after_MF = filter(g_arr, [1], recv_signal);
recv_noise_after_MF = filter(g_arr, [1], noises);

subplot(2,1,2);
stem(recv_signal_after_MF);
stem(recv_signal_after_MF + recv_noise_after_MF);

disp(['Energy_syms=', num2str(sum(transmit_delta_sequence.^2))]);
disp(['Energy_transmitted_signal=', num2str(sum(transmit_signal.^2))]);
disp(['channel Es/n0(dB) = ', num2str(0.5*db(1/(2*sigma^2)))]);

recv_delay = 16;
recv_sample_signal_after_MF = recv_signal_after_MF(recv_delay+1:oversample_rate:end);
recv_sample_noise_after_MF = recv_noise_after_MF(recv_delay+1:oversample_rate:end);
disp(['Energy_sampled_noise_after_MF=', num2str(sum(recv_sample_noise_after_MF.^2))]);
disp(['Energy_sampled_signal_after_MF=', num2str(sum(recv_sample_signal_after_MF.^2))]);
disp(['statistical SNR(dB) = ', num2str(...
    1/2*db(sum(recv_sample_signal_after_MF.^2)/sum(recv_sample_noise_after_MF.^2)))]);
% Expected: SNR <= Es/n0 + 3dB