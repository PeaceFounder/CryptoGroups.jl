# # Key Encpsulation Mechanism
# Essentailly Diffie-Hellman
using Test
using CryptoGroups
using Random: RandomDevice

function keygen(g::Group)

    sk = rand(RandomDevice(), 2:order(g) - 1)
    pk = octet(g^sk)

    return sk, pk
end

function encap(g::G, pk::Vector{UInt8}) where G <: Group

    y = G(pk)

    r = rand(RandomDevice(), 2:order(G) - 1) 
    t = y^r
    k = octet(t) # further hashing recomended
    c = g^r

    return k, c
end

function decap(sk, c::Group)

    t = c^sk
    k = octet(t)

    return k
end


g = @ECGroup{P_192}()
sk, pk = keygen(g)

k, c = encap(g, pk)

k′ = decap(sk, c)

@test k′ == k
