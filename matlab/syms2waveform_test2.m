%包含syms2waveform.m与waveform2syms.m的测试
clear all;
setup_wave;

%% Generate complex symbols with energy 1.
N_syms = 10;
N_sim = 1;

n0 = 2*0.1*0.1;
len_signal = N_syms * waveform_conf.oversample_rate + 2*waveform_conf.Group_delay;

for sim_iter = 1:N_sim
    syms_I = 1-2*(rand([N_syms, 1])>0.5);
    syms_Q = 1-2*(rand([N_syms, 1])>0.5);
    syms=syms_I.'+j*syms_Q.';
    [transmit_signal] = syms2waveform(syms,waveform_conf);
    
    [recv_syms] = waveform2syms(transmit_signal,n0,N_syms,waveform_conf);


end

%% Analyze RF signal spectrum.
% Spectrum estimation.
% L_spec = length(transmit_signal_RF{1});
% psd = zeros(L_spec,1);      % Calculate power spectrum density estimator.
% 
% for sim_iter = 1:N_sim
%     psd = psd + abs(fft(transmit_signal_RF{sim_iter})).^2;
% end
% psd = fftshift(psd / N_sim);
% figure(5);
% plot(psd);
