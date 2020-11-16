function [Sx, freq] = psd_est(x_mat, signal_len, fs, plot_switch)
% input matrix: 
% each row: 1 realization
% Assume: row(x_mat) = N_sampleLen >= signal_len.
% Using Period-Graph method.
    if ~exist('plot_switch', 'var') || isempty(plot_switch)
        plot_switch = false;
    end
    
    [M, N_sampleLen] = size(x_mat);
    Sx = zeros(N_sampleLen, 1);
    for k = 1:M
        Sx = Sx + (abs(fft(x_mat(k,:).'))).^2/signal_len;
    end
    Sx = Sx / M;    % calculate mean, reducing variance.
    
    % Convert digital PSD into analog PSD.
    PSD_len = ceil(signal_len/2);
    Sx = (Sx(1:PSD_len) / fs);
    freq = (0:PSD_len-1) / N_sampleLen * fs;
    
    if plot_switch
        plot(freq, Sx);
        title('Power Spectrum Density');
        xlabel('Frequency (Hz)');
        ylabel('PSD');
    end
    
end
