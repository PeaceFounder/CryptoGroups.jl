using ..Fields: Field
using ..CryptoGroups.Utils: static
using ..CryptoGroups: isstrict, order

abstract type EllipticCurve end
abstract type AbstractPoint end

(::Type{P})(x) where P <: AbstractPoint = convert(P, x)

function (::Type{P})(x, y) where P <: AbstractPoint

    F = field(P)
    
    _x = F(x) # convert(F, x) could also be used
    _y = F(y)

    return P(_x, _y)
end

Base.:-(u::P, v::P) where P <: AbstractPoint = u + (-v)

Base.isless(x::P, y::P) where P <: AbstractPoint = gx(x) == gx(y) ? gx(x) < gx(y) : gy(x) < gy(y)

function validate(x::AbstractPoint, order::Integer, cofactor::Integer)

    oncurve(x) || throw(ArgumentError("Point is not in curve"))
    x * cofactor != zero(x) || throw(ArgumentError("Point is in cofactor subgroup"))

    return
end

struct ECPoint{P<:AbstractPoint, S} <: AbstractPoint # The same contract is satisfied thus a subtype
    p::P

    function ECPoint{P, S}(x::P; allow_zero=false, skip_validation=false) where {P <: AbstractPoint, S}

        if iszero(x) && !allow_zero
            msg = "Constructing an offcurve element zero. Use `allow_zero` to hide this warning."
            if isstrict()
                error(msg)
            else
                @warn msg
            end
        end

        if !skip_validation
            # This allows specializing for Montgomery curves where clamping is used
            validate(x, S.order, S.cofactor)
        end

        # A test with cofactor also here
        new{P, S}(x)
    end

    function ECPoint(p::P, order::Integer, cofactor::Integer; name=nothing) where P <: AbstractPoint
        EP = concretize_type(ECPoint{P}, order, cofactor; name)
        return EP(p)
    end

    Base.zero(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = new{P, S}(zero(P))

    ECPoint{P, S}(x::F, y::F) where {P <: AbstractPoint, S, F <: Field} = ECPoint{P, S}(P(x, y))
end

function concretize_type(::Type{ECPoint{P}}, order::Integer, cofactor::Integer; name=nothing) where P <: AbstractPoint
    svars = static(; order, cofactor, name)
    return ECPoint{P, svars}
end

order(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = convert(Integer, S.order) 
order(::P) where P <: ECPoint = order(P)

cofactor(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = convert(Integer, S.cofactor) 
cofactor(::P) where P <: ECPoint = cofactor(P)

name(::Type{ECPoint}) = nothing
name(::Type{ECPoint{P}}) where P <: AbstractPoint = nothing
name(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = isnothing(S.name) ? nothing : convert(Symbol, S.name)


eq(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = eq(P)
field(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = field(P)

Base.zero(::P) where P <: ECPoint = zero(P)

function Base.:+(x::P, y::P) where P <: ECPoint 
    if iszero(x)
        return y
    elseif iszero(y)
        return x
    elseif x.p == y.p
        return P(double(x.p); skip_validation=true)
    elseif x.p == -y.p
        return zero(P)
    else
        return P(x.p + y.p; skip_validation=true)
    end
end

Base.:-(x::P) where P <: ECPoint = P(-x.p; skip_validation=true)


Base.:*(x::P, n::Integer) where P <: ECPoint = P(x.p * mod(n, order(P)); skip_validation=true)
Base.:*(n::Integer, x::ECPoint) = x * n

Base.convert(::Type{ECPoint{P, S}}, x::NTuple{2}; allow_zero=false) where {P <: AbstractPoint, S} = ECPoint{P, S}(convert(P, x); allow_zero)
Base.convert(::Type{P}, x::P) where P <: ECPoint = x 

oncurve(p::ECPoint) = oncurve(p.p)

Base.:(==)(x::ECPoint, y::ECPoint) = x.p == y.p

gx(p::ECPoint) = gx(p.p)
gy(p::ECPoint) = gy(p.p)

value(p::ECPoint) = value(p.p)
