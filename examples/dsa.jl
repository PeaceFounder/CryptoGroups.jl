using Test
using CryptoGroups
using CryptoGroups.Utils: modinv
using Random: RandomDevice

# DSA example

struct DSA
    r::Integer
    s::Integer
end

function keygen(g::Group)

    sk = rand(RandomDevice(), 2:order(g) - 1)
    pk = octet(g^sk)

    return sk, pk
end

function sign(e::Integer, g::Group, sk::Integer)

    n = order(g)

    k = rand(RandomDevice(), 2:order(g) - 1)

    R = g^k

    r = R % n
    s = modinv(k, n) * (e + sk * r) % n 

    if r == 0 || s == 0
        return sign(e, g, sk)
    else
        return DSA(r, s) 
    end
end

function verify(e::Integer, P::G, pk::Vector{UInt8}, sig::DSA) where {G <: Group}

    n = order(G)

    (; r, s) = sig
    
    @assert 1 < r < n - 1
    @assert 1 < s < n - 1

    Q = G(pk)
    c = modinv(s, n)
    
    u₁ = e * c % n
    u₂ = r * c % n

    # Special cases must be handled by the developer
    if u₁ == 0
        W = Q ^ u₂
    elseif u₂ == 0
        W = P ^ u₁
    else
        W = P ^ u₁ * Q ^ u₂ 
    end

    ν = W % n

    return ν == r
end


g = @ECGroup{P_192}() 
sk, pk = keygen(g)


m = 42

sig = sign(m, g, sk)

@test verify(m, g, pk, sig)
