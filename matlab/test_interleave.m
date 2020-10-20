clear;
setup_mapper;
mapping_conf.interleave=true;
mapping_conf.out='hard';
mapping_conf.alg = 'zf'; 
N_b = 10031;
random_bits = (rand([1,N_b])>0.5);
mapped_syms = bit_mapping(random_bits, mapping_conf);
demapped_bits = bit_demapping(mapped_syms, N_b, mapping_conf, ...
    ones([1,length(mapped_syms)]), ch_conf, 0.0001);
find(xor(demapped_bits, random_bits))