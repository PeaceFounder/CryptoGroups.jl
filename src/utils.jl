using Paillier 

rngprime(rng::AbstractRNG,N) = Paillier.nbit_prime_of_size(rng,N)
rngint(rng::AbstractRNG,N) = Paillier.n_bit_random_number(rng,N)

function hexstrtoint(str)
    hexstr = lowercase(join(split(replace(str,"\n"=>"")," ")))
    return parse(BigInt,str,base=16)
end

macro hex_str(str)
    esc(hexstrtoint(str))
end
