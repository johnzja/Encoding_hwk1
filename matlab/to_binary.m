function ret=to_binary(ind, bin_length)
% returns: binary form of integer ind. LSB last, MSB first.
    ret = false([1, bin_length]);
    for k = 1:bin_length
        ret(bin_length+1-k) = mod(ind, 2);
        ind = floor(ind/2);
    end
end