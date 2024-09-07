# # Key Encpsulation Mechanism

# A Key Encapsulation Mechanims allows (KEM) is a widelly used techinque to assymetrically encrypt messages that only recepient with knowledge of it's secret key can decrypt. In it's essence it is a Diffie-Hellman key computation where the secret key is used to encrypt the message with a symmetric cypher. With CryptoGroups we can write the key computation in a group agnostic way by defining three functions:

# 1. `keygen`: Generates a private-public key pair;
# 2. `encap`: Encapsulates a shared secret key using the recipient's public key;
# 3. `decap`: Decapsulates the shared secret key using the recipient's private key.

# At the end of the code we demonstrate how the code can be used with P-192 curve and verify that the derived keys match. Notice that methods are defined in a group agnostic way. Together with symmetric primitives defined in [Nettle.jl](https://github.com/JuliaCrypto/Nettle.jl) this can be an excellent starting point for implementing some assymetric encryption specifications in Julia.

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

function decap(c::Group, sk::Integer)

    t = c^sk
    k = octet(t)

    return k
end


g = @ECGroup{P_192}()
sk, pk = keygen(g)

k, c = encap(g, pk)

k′ = decap(c, sk)

@test k′ == k
