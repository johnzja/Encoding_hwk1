function [recv_syms] = waveform2syms(transmit_signal,n0,N_syms,waveform_conf)
%Summary of this function goes here
%载波波形转电平符号的一次实现
%输入：transmit_signal         载波波形           
%      n0                      噪声n0 
%      N_syms                  复电平符号个数
%      waveform_conf           波形参数
%输出：recv_syms               接收端复电平符号序列
%   Detailed explanation goes here
    disp_flag = true;    
    oversample_rate = waveform_conf.oversample_rate;
    fs = waveform_conf.fs;     % sample rate
    fc = waveform_conf.fc;      % carrier freq = 1850Hz    
    %Group_delay=waveform_conf.Group_delay;
    g_arr=waveform_conf.g_arr;
    sigma=sqrt(n0/2);
    len_signal=length(transmit_signal);
    recv_syms=zeros(1,N_syms);
    noises = sigma * randn(size(transmit_signal));
    recv_signal_baseband = (2*transmit_signal).*exp(-1j*2*pi*((0:len_signal-1))*fc/fs);
    recv_noise_baseband = (2*noises).*exp(-1j*2*pi*((0:len_signal-1))*fc/fs);
    recv_signal_after_MF = filter(g_arr, [1], recv_signal_baseband);                                 % g_arr filter is an LPF.
    recv_noise_after_MF = filter(g_arr, [1], recv_noise_baseband);
    
    recv_wave_after_MF=recv_signal_after_MF+recv_noise_after_MF;
    for k = 1:N_syms
        recv_syms(k)=recv_wave_after_MF(1+(k-1)*oversample_rate);
    end
    
    if disp_flag
        figure;
        stem(real(recv_signal_after_MF));
        hold on;
        stem(real(recv_signal_after_MF + recv_noise_after_MF));
        legend('recv BB', 'recv BB with noise');
    end
%     disp(['Energy_syms=', num2str(sum(abs(transmit_delta_sequence).^2))]);
%     disp(['Energy_transmitted_signal=', num2str(sum(abs(transmit_signal).^2))]);
%     disp(['channel Es/n0(dB) = ', num2str(0.5*db(0.5/(2*sigma^2)))]);
% 
%     recv_delay = Group_delay*2;
%     recv_sample_signal_after_MF = recv_signal_after_MF(recv_delay+1:oversample_rate:end);       % complex
%     recv_sample_noise_after_MF = recv_noise_after_MF(recv_delay+1:oversample_rate:end);         % complex
%     disp(['Energy_sampled_noise_after_MF=', num2str(sum(abs(recv_sample_noise_after_MF).^2))]);
%     disp(['Energy_sampled_signal_after_MF=', num2str(sum(abs(recv_sample_signal_after_MF).^2))]);
%     disp(['statistical SNR(dB) = ', num2str(...
%         1/2*db(sum(abs(recv_sample_signal_after_MF).^2)/sum(abs(recv_sample_noise_after_MF).^2)))]);
%     % Expected: SNR <= Es/n0 + 3dB
    
end

