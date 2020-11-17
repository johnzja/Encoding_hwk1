% Raised Cosine designer.
oversample_rate = 8;
N_filter_len = (2*10+1);
N_single_sided_filter_len = (N_filter_len-1)/2;

alpha = 1/2;
fs = 16000;
BW = 3400 - 300;
omega_c = pi*BW/fs/(1+alpha);
L = (pi/omega_c);           % H(e^j0)

N_GSample = 8000;           % Bigger, better.
% Plot graph of G.
figure(1);
hold on;
omega_arr = 2*pi*(-N_GSample:1:N_GSample).'/(2*N_GSample);
G_arr = zeros(size(omega_arr));
for k=1:length(omega_arr)
    G_arr(k) = get_G(omega_arr(k), L, omega_c, alpha);
end
plot(omega_arr, G_arr);
set(gca, 'xtick', -pi:pi/2:pi);
set(gca,'XTickLabel',{'-\pi', '-\pi/2','0', '\pi/2', '\pi'});
title('|G(e^j^\omega)|^2');

%% Perform IDTFT, adding window.
% Formula: 1/(2*pi)??sqrt(G(??))d??
N_gSample = 32;
g_arr = zeros(N_gSample*2+1, 1);    % Central point @ g_arr[32].
Group_delay = N_gSample;

for k = 0:N_gSample
    to_be_integrated = sqrt(G_arr) .* exp(1j*omega_arr*k);
    g_arr(N_gSample+1+k) = real(trapz(omega_arr, to_be_integrated)/(2*pi));
    if k ~= 0
        g_arr(N_gSample-k+1) = g_arr(N_gSample+1+k);  % Symmetric.
    end
end
% figure(4);
% freqz(g_arr,[1]);

g_arr = g_arr .* kaiser(length(g_arr),1);          % Adding kaiser window.
g_arr = g_arr / sqrt(sum(g_arr.^2));    % Energy normalization.

disp(sum(g_arr.^2));
figure(2);
plot((-N_gSample:1:N_gSample).', g_arr);
title('g[n]');
% figure(3);
% freqz(g_arr,[1]);

save('data/rcf.mat', 'g_arr', 'alpha', 'Group_delay');

function G = get_G(omega, L, omega_c, alpha)
    % Assume -pi <= omega <= pi.
    omega = abs(omega);
    omega_z0 = (1-alpha) * omega_c;
    omega_z1 = (1+alpha) * omega_c;
    if omega <= omega_z0
        G = L;
    elseif omega <= omega_z1
        G = L/2*(1+cos(pi*(omega-omega_z0)/(2*alpha*omega_c)));
    else
        G = 0;
    end
end