clear;
load('data/rcf.mat');

setup_encoder;
setup_mapper;
setup_wave;

N_info_bits = 8000;                                                     % 1kB file to transmit.
Ebn0_arr = [20, 15, 12.5, 10, 8, 6, 5, 4, 3, 2.5, 2, 1.5, 1, 0.5, 0, -0.5, -1, -1.5];     % Eb/n0 in dB.

Es = 0.5;                                   % 1 sym in waveform channel: 0.5 energy.
Eb = Es/mapping_conf.bps * (conv_encoder_conf.n/conv_encoder_conf.k);   % Notice: if hard-decision, Eb is different.
n0_arr = Eb ./ (10.^(Ebn0_arr/10));         % n0 in linear scale.
N_n0s = length(n0_arr);

L_encoded = zeros(N_n0s, 1);
L_transmit_syms = zeros(N_n0s, 1);
plot_iter = 0;

for n0_iter = 1:N_n0s
    random_bits = (rand([1, N_info_bits])>0.5);
    encoded_bits = conv_encode(random_bits, conv_encoder_conf);
    syms_transmit = bit_mapping(encoded_bits, mapping_conf);  
    wave_transmit = syms2waveform(syms_transmit, waveform_conf, false);
    [syms_recv, wave_recv] = waveform2syms(wave_transmit, n0_arr(n0_iter), length(syms_transmit), ...
            waveform_conf, false);
        
    % Draw Constellation for syms_recv.
    if plot_iter == 0
        figure;
    end
    
    subplot(2,2,plot_iter+1);
    axis equal;
    hold on;
    scatter(real(syms_recv), imag(syms_recv));
    scatter(real(syms_transmit), imag(syms_transmit));
    title(['Eb/n_0: ', num2str(Ebn0_arr(n0_iter)), 'dB']);
    legend('recv', 'transmit');
    
    if plot_iter == 3
        pause;
        close all;
    end
    plot_iter = mod(plot_iter+1, 4);
    
end


