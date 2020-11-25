encrypter.enable = true;
encrypter.method = 'DES';
if strcmp(encrypter.method, 'DES')
    encrypter.key = keygen();
    encrypter.blk_size_info_bits = 64;
    encrypter.blk_size_encrypted_bits = 64;
elseif strcmp(encrypter.method, 'RSA')
    encrypter.key = [988027, 283, 449467]; % [n,e,d]
    encrypter.blk_size_info_bits = floor(log2(encrypter.key(1)));
    encrypter.blk_size_encrypted_bits = ceil(log2(encrypter.key(1)));
end
