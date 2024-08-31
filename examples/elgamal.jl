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

# A scalar example with 3 row elements

g = @ECGroup{P_192}()

x = 42
pk = g^x

m = (g^2, g^3, g^4)

cyphertexts = enc(g, pk, m, (4, 5, 6))

plaintexts = dec(cyphertexts, x)

@test plaintexts == m #

# A vector case with 3 rows

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
