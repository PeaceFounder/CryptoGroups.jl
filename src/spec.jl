using .Utils: StaticBitVector, StaticBigInt, static, @check
using .Fields: F2GNB, F2PB, @F2PB, FP, Field, PrimeField, BinaryField, tobits, value, reducer
using .Curves: AbstractPoint, ECPoint, AffinePoint, Weierstrass, BinaryCurve, gx, gy, field, eq, a, b, cofactor
using .Specs: MODP, Koblitz, ECP, EC2N, GroupSpec, PB, GNB, curve

concretize_type(::Type{FP}, p::Integer) = FP{static(p)} 

concretize_type(::Type{BinaryCurve}, a::BitVector, b::BitVector) = BinaryCurve{StaticBitVector(a), StaticBitVector(b)}
concretize_type(::Type{BinaryCurve}, a::F, b::F) where F <: BinaryField = concretize_type(BinaryCurve, tobits(a), tobits(b)) 


static_int(x::Integer) = x
static_int(x::BigInt) = StaticBigInt(x)

concretize_type(::Type{Weierstrass}, a::Integer, b::Integer) = Weierstrass{static_int(a), static_int(b)}

concretize_type(::Type{Weierstrass}, a::BitVector, b::BitVector) = Weierstrass{StaticBitVector(a), StaticBitVector(b)}
concretize_type(::Type{Weierstrass}, a::F, b::F) where F <: BinaryField = concretize_type(Weierstrass, tobits(a), tobits(b))

concretize_type(::Type{P}, curve::Koblitz) where P <: AbstractPoint = concretize_type(P, curve.bec)
concretize_type(::Type{ECPoint{P}}, curve::Koblitz; name = name(curve)) where P <: AbstractPoint = concretize_type(ECPoint{P}, curve.bec; name)

function concretize_type(::Type{ECPoint{P}}, curve::GroupSpec; name = name(curve)) where P <: AbstractPoint

    Q = concretize_type(P, curve)

    _order = order(curve)
    _cofactor = cofactor(curve)
    
    R = concretize_type(ECPoint{Q}, _order, _cofactor; name)

    return R
end

concretize_type(::Type{ECPoint}, spec::GroupSpec; name = name(spec)) = concretize_type(ECPoint{AffinePoint}, spec; name)

concretize_type(::Type{F2GNB}, N::Int) = concretize_type(F2GNB, GNB(N))

function concretize_type(::Type{AffinePoint{Weierstrass, F}}, curve::ECP) where F <: PrimeField

    (; p, a, b) = curve

    P = AffinePoint{concretize_type(Weierstrass, a, b), concretize_type(F, p)}
    
    return P
end

concretize_type(::Type{AffinePoint}, spec::ECP) = concretize_type(AffinePoint{Weierstrass, FP}, spec)

function concretize_type(::Type{<:F2PB}, basis::PB)
    (; f) = basis
    f_bits = reverse(BitVector(i in f for i in 0:maximum(f)))
    return @F2PB{f_bits}
end


function concretize_type(::Type{<:F2GNB}, basis::GNB) #where F <: BinaryField
    (; m, T) = basis
    return F2GNB{m, T}
end

concretize_type(::Type{<:FP}, basis::Integer) = FP{static(basis)}


concretize_type(::Type{BinaryCurve}, curve::EC2N) = concretize_type(BinaryCurve, a(curve), b(curve))

function concretize_type(::Type{AffinePoint{BinaryCurve, F}}, curve) where F <: BinaryField
    P = AffinePoint{concretize_type(BinaryCurve, curve), concretize_type(F, curve.basis)}
    return P
end

concretize_type(::Type{AffinePoint}, curve::EC2N{GNB}) = concretize_type(AffinePoint{BinaryCurve, F2GNB}, curve)
concretize_type(::Type{AffinePoint}, curve::EC2N{PB}) = concretize_type(AffinePoint{BinaryCurve, F2PB}, curve)

concretize_type(::Type{AffinePoint{BinaryCurve, F}}, curve::Koblitz) where F <: BinaryField = concretize_type(AffinePoint, curve.bec)

function spec(::Type{P}; kwargs...) where P <: AbstractPoint

    n = hasmethod(order, Tuple{Type{P}}) ? order(P) : nothing
    h = hasmethod(cofactor, Tuple{Type{P}}) ? cofactor(P) : nothing

    EQ = eq(P)
    F = field(P)

    return spec(EQ, F; n, cofactor = h, kwargs...)
end

spec(p::P) where P <: AbstractPoint = spec(P; Gx=value(gx(p)), Gy=value(gy(p)))


function spec(::Type{EQ}, ::Type{F}; n=nothing, cofactor=nothing, Gx=nothing, Gy=nothing) where {EQ <: Weierstrass, F <: PrimeField}
    
    _a = a(EQ)
    _b = b(EQ)

    p = modulus(F)

    return ECP(p, n, _a, _b, cofactor, Gx, Gy) # I will need to add H in the end here
end


spec(::Type{F}) where F <: F2PB = PB(reducer(F))
spec(::Type{F2GNB{N, T}}) where {N, T} = GNB(N, T)


function spec(::Type{EQ}, ::Type{F}; n=nothing, cofactor=nothing, Gx=nothing, Gy=nothing) where {EQ <: BinaryCurve, F <: BinaryField}
    
    _a = convert(BitVector, a(EQ))
    _b = convert(BitVector, b(EQ))  ### Perhaps I would be better to do conversion at acessor methods!

    basis = spec(F)

    if isnothing(Gx) && isnothing(Gy)
        return EC2N(basis, n, _a, _b, cofactor) # I will need to add H in the end here
    else
        return EC2N(basis, n, _a, _b, cofactor, Gx, Gy) 
    end
end


concretize_type(::Type{ECGroup{P}}, spec::GroupSpec; name = name(spec)) where P <: ECPoint = ECGroup{concretize_type(P, spec; name)}
concretize_type(::Type{ECGroup}, spec::GroupSpec; name = name(spec)) = ECGroup{concretize_type(ECPoint, spec; name)}

concretize_type(::Type{PGroup}, p::Integer, q::Union{Integer, Nothing}; name::Union{Symbol, Nothing} = nothing) = PGroup{static(; p, q, name)}
concretize_type(::Type{PGroup}, spec::MODP; name = nothing) = concretize_type(PGroup, spec.p, spec.q; name)

"""
    spec(name::Symbol)::GroupSpec

Gets a group specification from a provided canonical `name` of the group specification which are hardcoded into the library. 

# Examples

```julia
# NIST elliptic curve P-192 specification:
p192_spec = spec(:P_192)
P192 = concretize_type(ECGroup, p192_spec)

# Modular prime group specification
modp_spec = spec(:RFC5114_2048_224)
MODP = concretize_type(PGroup, modp_spec)
```

See also `concretize_type`, `@ECGroup`, `@PGroup`
"""
function spec(x::Symbol)
    try
        return curve(x)
    catch
        return modp_spec(x)
    end
end


spec(::Type{ECGroup{P}}) where P = spec(P)
spec(g::ECGroup) = spec(g.x)

spec(::Type{G}) where G <: PGroup = MODP(; p = modulus(G), q = order(G))

(::Type{P})() where P <: ECPoint = P(generator(curve(name(P))))
(::Type{ECGroup{P}})() where P <: ECPoint = ECGroup{P}(P())
(::Type{G})() where G <: PGroup = G(generator(modp_spec(name(G))))


# Shall be added to CryptoGroups
iscompatable(x::GroupSpec, y::GroupSpec) = false

function iscompatable(x::ECP, y::ECP)

    x.p == y.p || return false
    x.n == y.n || return false
    x.a == y.a || return false
    x.b == y.b || return false
    x.cofactor == y.cofactor || return false

    return true
end

function iscompatable(x::EC2N, y::EC2N)

    x.basis == y.basis || return false
    x.n == y.n || return false
    x.a == y.a || return false
    x.b == y.b || return false
    x.cofactor == y.cofactor || return false

    return true
end

iscompatable(G::Type{<:Group}, Q::Type{<:Group}) = iscompatable(spec(G), spec(Q))


function _convert(G::Type{<:Group}, q::Group)
    @check iscompatable(G, typeof(q)) "Conversion not possible as group types represent different specifications"
    return G(octet(q))
end

Base.convert(G::Type{<:Group}, q::Group) = _convert(G, q)
Base.convert(G::Type{ECGroup{P}}, q::Group) where P <: ECPoint = _convert(G, q)
