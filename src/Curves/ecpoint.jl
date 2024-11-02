using ..Fields: Field
using ..CryptoGroups.Utils: static
using ..CryptoGroups: isstrict, order

abstract type EllipticCurve end

"""
    abstract type AbstractPoint end

An elliptic curve interface type. It is **not** expected that constructor of the group does instance validation and hence is unsafe. For safe use see `ECPoint`. The elliptic curve interface expects `*`, `+`, `-`, `gx`, `gy` and `zero` element to be implemented wheras the `order` and `cofactor` methods are added with the `ECPoint` subtype. In addition `oncurve` method is expected when used in combination with `ECPoint`.

It is expected that a `point::AbstractPoint` can be constructed back from it's coordinates and also from their octet representations:
```julia
P = typeof(point)
x, y = gx(point), gy(point)
point == P(x, y) == P(octet(x), octet(y))
```
Point compression is automatically supported through using `octet` on `gx(point)` and `gy(point)`. For point decompression from a compressed octet the user must implement a constructor for the point type `P(x::Vector{UInt8}, ỹ::Bool)` or a direct method `P(po::Vector{UInt8})`. 
"""
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

function validate(x::P, order::Integer, cofactor::Integer) where P <: AbstractPoint

    oncurve(x) || throw(ArgumentError("Point is not in curve"))
    #x * cofactor != zero(P) || throw(ArgumentError("Point is in cofactor subgroup"))
    if cofactor != 1
        !iszero(x * cofactor) || throw(ArgumentError("Point is in cofactor subgroup"))
    end

    return
end

"""
    rem(p::AbstractPoint, q::T)::T where T <: Integer

Computes a remainder of point `x` coordinate for a modulus `q`. The output is such that DSA and ECDSA of FIPS 186-4 standart
can be combined into single imlementation

"""
Base.rem(p::AbstractPoint, q::Integer) = rem(gx(p), q)


# TODO: At the moment multiplication with `0` is allowed which may be restricted in a future refactor. 

"""
    struct ECPoint{P<:AbstractPoint, S} <: AbstractPoint
        p::P
    end

A wrapper type for elliptic curve points ensuring their safe use. The constructor validates that the point is on the curve
and that it is not in cofactor subgroup. This check can be bypassed by specializing on `validate` method. Also adds safety 
checks when doing a point summation in case when points are equal, inverse of each other and treats zero properly. 

In addition to `AbstractPoint` methods `*`, `+`, `-`, `gx`, `gy` and `zero`, `ECPoint` implements `order`, `cofactor`
"""
struct ECPoint{P<:AbstractPoint, S} <: AbstractPoint # The same contract is satisfied thus a subtype
    p::P

    function ECPoint{P, S}(x::P; allow_zero=false, skip_validation=false) where {P <: AbstractPoint, S}

        if !skip_validation
            if iszero(x) 
                if !allow_zero
                    msg = "Constructing an offcurve element zero. Use `allow_zero` to hide this warning/error."
                    isstrict() ? throw(ArgumentError(msg)) : @warn msg
                end
            else 
                # This allows specializing for Montgomery curves where clamping is used
                validate(x, S.order, S.cofactor)
            end
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

ECPoint{P, S}(x; allow_zero=false, skip_validation=false) where {P <: AbstractPoint, S} = ECPoint{P, S}(convert(P, x); allow_zero, skip_validation)

Base.iszero(x::ECPoint) = iszero(x.p)

"""
    concretize_type(::Type{ECPoint{P}}, order::Integer, cofactor::Integer; name=nothing) where P <: AbstractPoint

Constructs an `ECPoint` for a concrete point type `P <: AbstractPoint` with specified order and cofactor. The name is
used for aliasing group parameters with `@ECPoint{$name}` display semantics. 
"""
function concretize_type(::Type{ECPoint{P}}, order::Integer, cofactor::Integer; name=nothing) where P <: AbstractPoint
    svars = static(; order, cofactor, name)
    return ECPoint{P, svars}
end

"""
    order(::Union{P, Type{P}}) where P <: ECPoint

Gets an order of a point type or instance. 
"""
order(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = convert(Integer, S.order) 
order(::P) where P <: ECPoint = order(P)

"""
    cofactor(::Union{P, Type{P}}) where P <: ECPoint

Gets a cofactor of a point type or instance. 
"""
cofactor(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = convert(Integer, S.cofactor) 
cofactor(::P) where P <: ECPoint = cofactor(P)

name(::Type{ECPoint}) = nothing
name(::Type{ECPoint{P}}) where P <: AbstractPoint = nothing
name(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = isnothing(S.name) ? nothing : convert(Symbol, S.name)

eq(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = eq(P)
field(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = field(P)
field(::Type{ECPoint{P}}) where {P <: AbstractPoint} = field(P)

"""
    zero(::Union{P, Type{P}}) where P <: AbstractPoint

Constructs a zero element representation of a point type. Implementation of arithmetic can be delegated to `ECPoint` wrapper
that treats the cases safely.
"""
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

function Base.:*(x::P, n::Integer; skip_validation = false) where P <: ECPoint 

    n_mod = mod(n, order(P))

    if !skip_validation && n_mod == 0
        msg = "A point is being multiplied with `0`" 
        if isstrict()
            error(msg)
        else
            @warn msg
            return zero(P)
        end
    end

    return P(x.p * n_mod; skip_validation=true)
end

Base.:*(n::Integer, x::ECPoint; skip_validation = false) = *(x, n; skip_validation)

Base.convert(::Type{ECPoint{P, S}}, x; allow_zero=false) where {P <: AbstractPoint, S} = ECPoint{P, S}(convert(P, x); allow_zero)
Base.convert(::Type{P}, x::P) where P <: ECPoint = x 

"""
    oncurve(p::AbstractPoint)::Bool

Checks if a point coordinates satisfy elliptic curve equation.
"""
oncurve(p::ECPoint) = oncurve(p.p)

Base.:(==)(x::ECPoint, y::ECPoint) = x.p == y.p

"""
    gx(p::AbstractPoint)::Field

Returns `x` coordinate of a point instance. See also `gy`, `value`, `octet`.
"""
gx(p::ECPoint) = gx(p.p)

"""
    gy(p::AbstractPoint)::Field

Returns `y` coordinate of a point instance. See also `gx`, `value`, `octet`.
"""
gy(p::ECPoint) = gy(p.p)

value(p::ECPoint) = value(p.p)
