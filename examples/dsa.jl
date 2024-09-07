# # Digital Signature Algorithm

# This Julia code demonstrates a group-agnostic implementation of the Digital Signature Algorithm (DSA), a widely used standard for message authentication. The implementation defines three key functions: keygen for key generation, sign for creating signatures, and verify for signature verification. While it omits the message hashing step for simplicity, the code adheres to FIPS 186-4 compliance, with safety checks integrated into the group element constructors. 

# The implementation is designed to work with any abstract group, showcasing Julia's multiple dispatch capabilities and allowing for flexibility across different cryptographic settings. To illustrate its functionality, the code includes a test case using the P-192 elliptic curve group. This example serves as a foundational demonstration of DSA and is foundation of [CryptoSignatures.jl](https://github.com/PeaceFounder/CryptoSignatures.jl) which offers extensive test suite and uses deterministic randomness for $k$ generation.

using Test
using CryptoGroups
using CryptoGroups.Utils: @check
using Random: RandomDevice

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
    s = invmod(k, n) * (e + sk * r) % n 

    if r == 0 || s == 0
        return sign(e, g, sk)
    else
        return DSA(r, s) 
    end
end

function verify(e::Integer, P::G, pk::Vector{UInt8}, sig::DSA) where {G <: Group}

    n = order(G)

    (; r, s) = sig
    
    @check 1 < r < n - 1
    @check 1 < s < n - 1

    Q = G(pk)
    c = invmod(s, n)
    
    u₁ = e * c % n
    u₂ = r * c % n

    ## Special cases must be handled by the developer
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
