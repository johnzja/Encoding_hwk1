function ret = from_binary(bin_form)
% returns: integer of which the binary form is bin_form.
% Notice: LSB last, MSB first.
    L=length(bin_form);
    nv = (2.^(L-1:-1:0)).';
    ret = bin_form*nv;
end