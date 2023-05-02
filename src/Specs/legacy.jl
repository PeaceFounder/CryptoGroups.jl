### Some prime group generation algorithms. References:
# + https://crypto.stackexchange.com/questions/820/how-does-one-calculate-a-primitive-root-for-diffie-hellman
# + https://math.stackexchange.com/questions/124408/finding-a-primitive-root-of-a-prime-number

# Eventually will need to revisit this implmentation to comply with CryptoGroups standart

using Primes: isprime, nextprime
using Random: AbstractRNG, default_rng

function n_bit_random_number(rng::AbstractRNG, len::Integer)
    max_n = ( BigInt(1) << len ) - 1
    if len > 2
        min_n = BigInt(1) << (len - 1)
        return rand(rng, min_n:max_n)
    end
    return rand(rng, 1:max_n)
end

function nbit_prime_of_size(rng::AbstractRNG, n_bits::Integer)
    # generate a random nbit number
    r = n_bit_random_number(rng, n_bits)
    return nextprime(r)
end


rngprime(rng::AbstractRNG, N) = nbit_prime_of_size(rng, N)
rngint(rng::AbstractRNG, N) = n_bit_random_number(rng, N)

### I could manage to get 512 security
# Something wrong with this one
function sophie_germain_group(rng::AbstractRNG, g::Integer, t::Int)
    @assert g!=1
    while true
        q = rngprime(rng, 2*t)
        p = 2*q + 1
        if g!=p-1 && isprime(p)
            return MODP(; p, g, q)
        end
    end
end

function generate_qp(rng::AbstractRNG, t)
    while true
        q = rngprime(rng, t)
        #q == 1 && continue

        r = rngint(rng, t)

        p = q*r + 1
        if isprime(p)
            return q, p
        end
    end
end

generate_qp(t::Int) = generate_qp(default_rng(), t)


function dsa_standart_group(rng::AbstractRNG, tp::Int, tg::Int)
    q, p = generate_qp(rng, 2*tp)
    while true
        u = mod(rngint(rng, tg), p)
        g = powermod(u, div((p-1), q), p)

        if g!=1 && g!=0
            return MODP(; p, g, q)
        end
    end
end

dsa_standart_group(tp::Int, tg::Int) = dsa_standart_group(default_rng(), tp, tg)


function generate_g(p::Integer, q::Integer; seed = Vector{UInt8}("SEED"), nr = 10, hasher = "sha256")

    prg = PRG(hasher; s = seed)
    sp = MODP(;p, q)

    h = rand(prg, sp, 1; nr)

    rand(prg, sp, 10; nr)

    return h[1]
end
