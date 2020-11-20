function [K] = keygen()
    K = randi(2,[8,8]) - 1;
    for row=1:8
        K(row,8) = mod(sum(K(row,1:7)),2);
    end
    K = reshape(K',[1,64]);
end

