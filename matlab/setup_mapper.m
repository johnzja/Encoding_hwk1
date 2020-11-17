%% Config all the parameters.
ch_conf.b = 0.5;
ch_conf.rho = 0.95;
mapping_conf.M = [1+1j, 1-1j, -1+1j, -1-1j]./sqrt(2);  % grey code with Ps=1.
                                    % [-1,1] [1, 1j, -1j, -1]
                                    % [1, 1/sqrt(2)+1j/sqrt(2), -1/sqrt(2)+1j/sqrt(2), ...
                                    % j, 1/sqrt(2)-1j/sqrt(2), -i, -1, -1/sqrt(2)-1j/sqrt(2) ]
mapping_conf.bps = log2(length(mapping_conf.M));
mapping_conf.mode = 'ch_known';     % ['ch_unknown', 'ch_known']

mapping_conf.alg = 'zf';            % ['zf', 'mse'];
mapping_conf.out = 'L2';            % ['hard', 'LLR', 'L2']. This parameter is automatically set in joint-simulation.
mapping_conf.pilotrate = 4;         % one pilot every N symbols
mapping_conf.est = 'kalman';        % ['linear','quad','kalman']. Methods for estimating channels
mapping_conf.interleave=false;      % interleave
mapping_conf.depth=8;               % interleaving depth