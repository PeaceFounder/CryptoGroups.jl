using CryptoGroups

struct ElGamalRow{G<:Group, N}
    a::NTuple{N, G}
    b::NTuple{N, G}
end

ElGamalRow(m::Vector{NTuple{N, G}}) where {N, G<:Group} = ElGamalRow(one.(m), m)

enc(g::G, pk::G, e::ElGamalRow{G, N}, r::NTuple{N, <:Integer}) where {N, G <: Group} = ElGamalRow(e.a .* g .^ r, e.b .* pk .^ r)
enc(g::G, pk::G, m::NTuple{N, G}, r::NTuple{N, <:Integer}) where {N, G <: Group} = enc(g, pk, ElGamalRow(m), r)

dec(e::ElGamalRow{G}, x::Integer) = e.b ./ e.a .^ x


g = @ECGroup{P_192}()

x = 42
pk = g^x

m = (g^2, g^3, g^4)

cyphertexts = enc(g, pk, m, (4, 5, 6))

plaintexts = dec(cyphertexts, x)

@test plaintexts == m #
