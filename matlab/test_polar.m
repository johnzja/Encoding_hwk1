addpath('utils/');

design_epsilon = 0.32;
crc_length = 16;
[gen, det, g] = get_crc_objective(crc_length);
n = 10;
N = 2^n;
K = N/2 + crc_length;   % Final coding rate = 1/2.
ebno_vec = [2 2.5];     % row vec, you can write it like [1 1.5 2 2.5 3] 
list_vec = [1 16];      % row vec, you can write it like [1 4 16 32 ...]. The first element is always 1 for acceleration purpose. The ramaining elements are power of two.
max_runs = 1e7;
max_err = 100;
resolution = 1e4;       % the results are shown per max_runs/resolution.
[bler, ber] = simulation(N, K, design_epsilon, max_runs, max_err, resolution,  ebno_vec, list_vec, gen, det, g, crc_length);