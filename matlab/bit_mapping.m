function syms=bit_mapping(bit_stream, mapping_conf)
% mapping_conf(struct):
    mapping_vector = mapping_conf.M;
    bit_per_symbol = mapping_conf.bps;
    mode = mapping_conf.mode;
    
    L = length(bit_stream);
    r = mod(L, bit_per_symbol);
    if r
        bit_stream = [bit_stream, false([1, bit_per_symbol-r])];%补零
        L = length(bit_stream);
    end
    Nb = L; Ns = Nb/bit_per_symbol;%计算需要的符号个数Ns
    syms = zeros([1, Ns]);
    
    nv = (2.^(bit_per_symbol-1:-1:0)).';
    for k=1:Ns
        index = 1+(k-1)*bit_per_symbol;
        slice = bit_stream(index:index+bit_per_symbol-1);
        n = slice*nv+1;
        syms(k) = mapping_vector(n);
    end
    %interleave
    interleave=mapping_conf.interleave;
    depth=mapping_conf.depth;
    if interleave
        Ls=length(syms);
        block=ceil(Ls/depth);
        remainder=mod(Ls, depth);
        syms_interleave=zeros(size(syms));
        for D=1:depth
            for B=1:block
                if(D<=remainder)
                    syms_interleave((D-1)*block+B)=syms((B-1)*depth+D);
                elseif((B-1)*depth+D<=Ls)
                    syms_interleave((D-1)*block+B-(D-remainder-1))=syms((B-1)*depth+D);
                end
            end
        end
        syms=syms_interleave;
    end
    
    % add pilots if we assume no CSI at Rx
    if strcmp(mode, 'ch_unknown')
        pilot_rate = mapping_conf.pilotrate;
        syms_with_pilots = zeros([1, Ns + ceil(Ns/pilot_rate)]);
        pnt = 1;
        for k=1:Ns
            if mod(pnt,pilot_rate) == 1
                syms_with_pilots(pnt) = 1.0;
                pnt = pnt + 1;
            end
            syms_with_pilots(pnt) = syms(k);
            pnt = pnt + 1;
        end
        syms_with_pilots = [syms_with_pilots 1.0];
        syms = syms_with_pilots;
    end
end