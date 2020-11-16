load('data/rcf.mat');
waveform_conf.oversample_rate = 8;
waveform_conf.fs = 16000;     % sample rate
waveform_conf.fc = 1850;      % carrier freq = 1850Hz
waveform_conf.g_arr=g_arr;
waveform_conf.Group_delay=Group_delay;
