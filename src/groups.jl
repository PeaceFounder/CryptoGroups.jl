using .Utils: int2octet, octet2int
import .Fields: value, modulus, octet
using .Curves: ECPoint, gx, gy

"""
    abstract type Group end

A cyclic group interface type. It is expected that constructor of the group does instance validation and throws an `ArgumentError`
if it can't be constructed. Identity element can be constructed with explicit argument `allow_one=true` and validation can be opted out with `skip_validation=true` argument for the constructor which can be useful for deserialization from a secure location locally. 

The group element shall support construction from vector of bytes (octet) as well as `octet` method itself. In addition `value` method is expected that provides most simple representation useful for debugging purposes. Given `g::Group` a following constructor methods are supported:
```julia
G = typeof(g)
g == G(octet(g)) == G(value(g))
one(G) == G(octet(one(G)), allow_one=true)
```

The group must support `order`, `*`, `^`, `rem` and identity element via `one`. For concrete implementations see `PGroup` and `ECGroup`. 
"""
abstract type Group end

Base.broadcasted(f::Function, x::Group, y::AbstractVector{<:Integer}) = f.((x for i in 1:length(y)), y)
Base.broadcasted(f::Function, x::Group, y::NTuple{N, <:Integer}) where N = ntuple(i -> f(x, y[i]), N)

Base.broadcasted(f::Function, x::G, y::NTuple{N, G}) where {N, G <: Group} = ntuple(i -> f(x, y[i]), N)
Base.broadcasted(f::Function, x::NTuple{N, G}, y::G) where {N, G <: Group} = Base.broadcasted(f, y, x)

Base.broadcasted(f::Function, x::G, y::Vector{G}) where G <: Group = f.((x for i in 1:length(y)), y)
Base.broadcasted(f::Function, x::Vector{G}, y::G) where G <: Group = Base.broadcasted(f, y, x)


"""
    value(g::Group)::Union{BigInt, Tuple{BigInt, BigInt}, Tuple{BitVector, BitVector}}

Converts the group instance to the most simple representation. Generally intended for debugging purposes and for serialization one should use `octet` as it also ensures consistent padding. 

See also `octet`, `Curves.gx`, `Curves.gy`
"""
function value(::Group) end

"""
    order(::Union{G, Type{G})::BigInt where G <: Group

Get the order of the group
"""
order(x::G) where G <: Group = order(G)

"""
    inv(g::G)::G where G <: Group

Computes inverse of the group element so that `inv(g) * g == one(G)`
"""
Base.inv(g::G) where G <: Group = g^(order(G) - 1) 

#import Base./

"""
    /(x::G, y::G)::G where G <: Group

Convinience method for `x * inv(y)`
"""
Base.:/(x::G, y::G) where G <: Group = x * inv(y)

name(x::G) where G <: Group = name(G)

"""
    iscompressable(g::Group)::Bool

Returns `true` if `octet` function accepts `mode=:compressed` as is groups based on eliptic curves. Currently for binary curves decompression is not implemented and hence this method returns false for them. 
"""
iscompressable(g::Group) = false 

"""
    struct ECGroup{P<:ECPoint} <: Group
        x::P
    end

Ellitpic curve group that maps `ECPoint` add and multiply to multiply and exponent semantics. The group can be instantiated from
`x`, `y` coordinates `ECGroup{P}(x, y)` which are passed to corepsonding elliptic curve point constructor or via octet representation as specified in FIPS 186-4 standart. To construct the type use `concretize_type` or `@ECGroup` macro.

# Example

```julia
p192_spec = spec(:P_192)
G = concretize_type(ECGroup, p192_spec)

# A macro syntax
G = @ECGroup{P_192}
```
If a group is instantiated from an existing specification then `G()` creates an instance of generator. A following ilustrates typical operations:

```julia
G = @ECGroup{P_192}
g = G()
g^3 * g^5/g^2 == (g^3)^2 == g^6
g^(order(G) - 1) * g == one(G)
one(G) * g == g
G(octet(g)) == G(value(g)) == g
```
See also `octet`, `value`, `concretize_type`, `spec`

"""
struct ECGroup{P<:ECPoint} <: Group
    x::P
end

ECGroup{P}(x, y) where {P <: ECPoint} = ECGroup{P}(P(x, y))
ECGroup{P}(x::Union{AbstractVector{UInt8}, NTuple{2, <:Union{Integer, BitVector}}}; allow_one=false, skip_validation=false) where {P <: ECPoint} = ECGroup{P}(P(x; allow_zero=allow_one, skip_validation))

Curves.field(::Type{ECGroup{P}}) where P <: ECPoint = field(P)

"""
    convert(::Type{G}, x; allow_one=false)::G where G <: Group

Converts representation of `x` into a group element type `G`. The conversion is safe as the validation checks are performed within
group constructors. In case identity element needs to be read into `allow_one` flag can be used to allow constructing identity elements. 
"""
Base.convert(::Type{ECGroup{P}}, x; allow_one=false) where P <: ECPoint = ECGroup{P}(convert(P, x; allow_zero=allow_one))
Base.convert(::Type{G}, x::G) where G <: Group = x


"""
    octet(g::ECGroup; mode = :uncompressed|:hybrid|:compressed)::Vector{UInt8}

Converts elliptic curve point to an octet representation as specified in FIPS 186-4 standart. The `mode=:uncompressed|:hybrid|:compressed` can be used to specify compression mode. 

See also `iscompressable`
"""
octet(g::ECGroup; mode = :uncompressed) = octet(g.x; mode)

"""
    *(x::G, y::G)::G where G <: Group

Multiplies two group elements.
"""
Base.:*(x::G, y::G) where G <: ECGroup = G(x.x + y.x)

"""
    ^(x::G, n::Integer)::G where G <: Group

Exponentiates the group element. In case `mod(n, order(G)) == 0` throws an error in strict mode or shows a warning in relaxed mode.

See also `isstrict` and `set_strict_mode`
"""
function Base.:^(x::G, n::Integer) where G <: ECGroup

    n_mod = mod(n, order(G))

    if n_mod == 0
        msg = "A bad exponent"
        if isstrict()
            error(msg)
        else
            @warn msg
            return one(G)
        end
    end
    
    if isone(x)
        return x
    else
        return G(*(n_mod, x.x; skip_validation = true))
    end
end

order(::Type{ECGroup{P}}) where P = order(P)

Base.isvalid(x::ECGroup) = isvalid(x.x)

Base.:(==)(x::G, y::G) where G <: ECGroup = x.x == y.x

modulus(::Type{ECGroup{P}}) where P <: ECPoint = modulus(P)

name(::Type{ECGroup}) = nothing
name(::Type{ECGroup{P}}) where P <: ECPoint = name(P)

Base.isless(x::G, y::G) where G <: ECGroup = isless(x.x, y.x)

Curves.gx(g::ECGroup) = gx(g.x)
Curves.gy(g::ECGroup) = gy(g.x)

"""
    one(::Type{G}) where G

Construct an identity element of the group.
"""
Base.one(g::ECGroup{P}) where P <: ECPoint = ECGroup{P}(zero(P))
Base.one(::Type{ECGroup{P}}) where P <: ECPoint = ECGroup{P}(zero(P))

value(g::ECGroup) = value(g.x)

iscompressable(g::ECGroup) = iscompressable(g.x)

"""
    rem(x::ECGroup, q::T)::T where T <: Integer

Computes remainder of the elliptic curve point
"""
Base.rem(x::ECGroup, q::Integer) = rem(x.x, q)

"""
    struct PGroup{S} <: Group
        g::BigInt
    end

Modulus prime group where `S` is a static type parameter encoding group properties. To instantiate a concrete type use a macor `@PGroup` or use `concretize_type`.

# Example

```julia
# Directly passing type arguments
G = concretize_type(PGroup, 23, 11) # where modulus is 23 and order is 11

# Using a macro for user specification
G = @PGroup{p = 23, q = 11}

# Using existing specification
modp = spec(:RFC5114_2048_224)
G = concretize_type(PGroup, modp)

# Using a macro for existing specification
G = @PGroup{RFC5114_2048_224}
```

If a group is instantiated from an existing specification then `G()` creates an instance of generator. A following ilustrates typical operations:

```julia
G = @PGroup{p = 23, q = 11}
g = G(2)
g^3 * g^5/g^2 == (g^3)^2 == g^6
g^(order(G) - 1) * g == one(G)
one(G) * g == g
G(octet(g)) == G(value(g)) == g
```
See also `octet`, `value`, `concretize_type`, `spec`
"""
struct PGroup{S} <: Group
    g::BigInt

    function PGroup{S}(x_int::Integer; allow_one::Bool=false, skip_validation=false) where S

        x = convert(BigInt, x_int)
        _order = order(PGroup{S})
        _modulus = modulus(PGroup{S})

        if !skip_validation
            if x == 1 
                if !allow_one
                    msg = "Constructing a degenerate element. Use `allow_one` to hide this warning"
                    isstrict() ? throw(ArgumentError(msg)) : @warn msg
                end
            elseif !isnothing(_order)
                0 < x < S.p || throw(ArgumentError("Element $x is not in range of the prime group with modulus $_modulus"))
                powermod(x, _order, _modulus) == 1 || throw(ArgumentError("Element $x is not an element of prime group with order $_order and modulus $_modulus"))
            end            
        end

        new{S}(x)
    end

    PGroup{S}(x::Vector{UInt8}; allow_one::Bool=false, skip_validation=false) where S = PGroup{S}(octet2int(x); allow_one, skip_validation)

    Base.one(::Type{PGroup{S}}) where S = new{S}(1)
end

Base.one(::PGroup{S}) where S = one(PGroup{S})


"""
    octet(x::PGroup)::Vector{UInt8}

Converts modulus prime group element into octet representation. A padding is added to match the length of modulus.
"""
octet(x::PGroup) = int2octet(value(x), bitlength(modulus(x)))

Base.convert(group::Type{<:PGroup}, element::Vector{UInt8}; allow_one::Bool=false) = convert(group, octet2int(element); allow_one)

"""
    modulus(::Union{G, Type{G}})::BigInt where G <: PGroup

Modulus of a prime group. It is not recommended to depend on this method in the codebase as it destroys polymorphism.
"""
modulus(::Type{PGroup{S}}) where S = BigInt(S.p)
modulus(::G) where G <: PGroup = modulus(G)

order(::Type{PGroup{S}}) where S = S.q isa Nothing ? nothing : BigInt(S.q)

name(::Type{PGroup}) = nothing

name(::Type{PGroup{S}}) where S = !(@isdefined S) || isnothing(S.name) ? nothing : convert(Symbol, S.name)


value(g::PGroup) = g.g

Base.convert(::Type{P}, x::Integer; allow_one=false) where P <: PGroup = P(BigInt(x); allow_one)

Base.isvalid(g::G) where G <: PGroup = value(g) != 1 && powermod(value(g), order(G), modulus(G)) == 1


import Base.*
function *(x::G, y::G) where G <: PGroup 
    if isone(x) 
        return y
    elseif isone(y)
        return x
    else
        return  G(mod(value(x) * value(y), modulus(G)); skip_validation=true)
    end
end


import Base.^
function ^(x::G, n::Integer) where G <: PGroup 

    n_mod = mod(n, order(G))

    # @assert n_mod != 0 "A bad exponent" 
    if n_mod == 0
        msg = "A bad exponent" 
        if isstrict()
            error(msg)
        else
            @warn msg
            return one(G)
        end
    end

    if isone(x)
        return x
    else
        return G(powermod(value(x), n_mod, modulus(G)); skip_validation=true)
    end
end


Base.inv(x::G) where G <: PGroup = G(invmod(value(x), modulus(G)); skip_validation=true)

Base.:(==)(x::G, y::G) where G <: PGroup = x.g == y.g

Base.isless(x::G, y::G) where G <: PGroup = value(x) < value(y)

"""
    rem(x::PGroup, q::T)::T where T <: Integer

Computes remainder of a prime group integer value
"""
Base.rem(x::PGroup, q::Integer) = rem(value(x), q)
