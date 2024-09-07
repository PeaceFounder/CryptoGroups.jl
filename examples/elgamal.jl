# # ElGamal Cryptosystem

# The ElGamal cryptosystem, as demonstrated in this Julia code example, is a fundamental component of many verifiable e-voting systems. This implementation showcases the core operations of the ElGamal encryption scheme, including key generation, encryption, and decryption and is implemented in a group agnostic way. 

# The example illustrates both scalar and vector cases of ElGamal encryption. It defines an `ElGamalRow` struct to represent encrypted messages, and implements functions for encryption (enc) and decryption (dec). The code demonstrates how to encrypt multiple messages simultaneously and how to work with tuples of group elements as well as vectors. By using Julia's multiple dispatch capabilities, the implementation provides a flexible and type-safe approach to handling different input types. The test cases at the end of the file are made using P-192 elliptic curve group verify the correctness of the encryption and decryption operations for both single-row and multi-row scenarios.

using Test
using CryptoGroups

struct ElGamalRow{G<:Group, N}
    a::NTuple{N, G}
    b::NTuple{N, G}
end

ElGamalRow(m::NTuple{N, G}) where {N, G<:Group} = convert(ElGamalRow{G, N}, m)
Base.convert(::Type{ElGamalRow{G, N}}, m::NTuple{N, G}) where {N, G <: Group} = ElGamalRow{G, N}(one.(m), m)

enc(g::G, pk::G, e::ElGamalRow{G, N}, r::NTuple{N, <:Integer}) where {N, G <: Group} = 
    ElGamalRow(e.a .* g .^ r, e.b .* pk .^ r)

enc(g::G, pk::G, m::NTuple{N, G}, r::NTuple{N, <:Integer}) where {N, G <: Group} = enc(g, pk, ElGamalRow(m), r)

dec(e::ElGamalRow{G}, x::Integer) where G <: Group = e.b ./ e.a .^ x

## A scalar example with 3 row elements

g = @ECGroup{P_192}()

x = 42
pk = g^x

m = (g^2, g^3, g^4)

cyphertexts = enc(g, pk, m, (4, 5, 6))

plaintexts = dec(cyphertexts, x)

@test plaintexts == m #

## A vector case with 3 rows

mvec = NTuple{3, typeof(g)}[
    (g^2, g^3, g^4),
    (g^4, g^7, g^8),
    (g^3, g^9, g^2)
]

rvec = NTuple{3, Int}[
    (4, 5, 6),
    (4, 5, 6),
    (4, 5, 6)
]

cyphertexts_vec = ((m, r) -> enc(g, pk, m, r)).(mvec, rvec)

plaintexts_vec = dec.(cyphertexts_vec, x)

@test plaintexts_vec == mvec
