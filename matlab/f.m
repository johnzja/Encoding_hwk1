function [output] = f(R, K, S)
    E = [32,1,2,3,4,5,4,5,6,7,8,9,8,9,10,11,12,13,...
        12,13,14,15,16,17,16,17,18,19,20,21,20,21,22,23,24,25,...
        24,25,26,27,28,29,28,29,30,31,32,1];
    P = [16,7,20,21,29,12,28,17,1,15,23,26,...
        5,18,31,10,2,8,24,14,...
        32,27,3,9,19,13,30,6,...
        22,11,4,25];
    R = mod(R(E)+K,2);
    R_S = zeros(1,32);
    for row = 1:8
        STable = S{row};
        blockOffset = 6*row-6;
        ind1 = 2*R(1,blockOffset + 1) + R(1,blockOffset + 5) + 1;
        ind2 = 8*R(1,blockOffset + 2) + 4*R(1,blockOffset + 3) +...
            2*R(1,blockOffset + 4) + R(1,blockOffset + 5) + 1;
        R_S(4*row-3: 4*row) = STable{ind1,ind2};
    end
    output = R_S(P);
end

