using .Utils: int2octet, octet2int
import .Fields: value, modulus, octet
using .Curves: ECPoint, gx, gy
using .Specs: PRG

abstract type Group end

Base.broadcasted(f::Function, x::Group, y::AbstractVector{<:Integer}) = f.((x for i in 1:length(y)), y)
Base.broadcasted(f::Function, x::Group, y::NTuple{N, <:Integer}) where N = ntuple(i -> f(x, y[i]), N)

Base.broadcasted(f::Function, x::G, y::NTuple{N, G}) where {N, G <: Group} = ntuple(i -> f(x, y[i]), N)
Base.broadcasted(f::Function, x::NTuple{N, G}, y::G) where {N, G <: Group} = Base.broadcasted(f, y, x)

Base.broadcasted(f::Function, x::G, y::Vector{G}) where G <: Group = f.((x for i in 1:length(y)), y)
Base.broadcasted(f::Function, x::Vector{G}, y::G) where G <: Group = Base.broadcasted(f, y, x)

#Base.broadcasted(::typeof(*), x::Vector{G}, y::G) where G <: Group =  x .* (y for i in 1:length(x))

order(x::G) where G <: Group = order(G)

Base.inv(g::G) where G <: Group = g^(order(G) - 1) 

import Base./
/(x::G, y::G) where G <: Group = x * inv(y)

name(x::G) where G <: Group = name(G)

Base.rand(prg::PRG, ::Type{G}, N::Integer; nr::Integer = 0) where G <: Group = convert(Vector{G}, rand(prg, spec(G), N; nr))

Base.ones(x::Vector{G}) where G <: Group = [one(i) for i in x] # need to extend this

iscompressable(g::Group) = false 

struct ECGroup{P<:ECPoint} <: Group
    x::P
end

ECGroup{P}(x, y) where {P <: ECPoint} = ECGroup{P}(P(x, y))
ECGroup{P}(x::Vector{UInt8}) where {P <: ECPoint} = ECGroup{P}(P(x))

Base.convert(::Type{ECGroup{P}}, x; allow_one=false) where P <: ECPoint = ECGroup{P}(convert(P, x; allow_zero=allow_one))
Base.convert(::Type{G}, x::G) where G <: ECGroup = x

octet(g::ECGroup; mode = :uncompressed) = octet(g.x; mode)

Base.:*(x::G, y::G) where G <: ECGroup = G(x.x + y.x)

function Base.:^(x::G, n::Integer) where G <: ECGroup

    n_mod = mod(n, order(G))

    # @assert n_mod != 0 "A bad exponent"    
    if n_mod == 0
        msg = "A bad exponent"
        if isstrict()
            error(msg)
        else
            @warn msg
        end
        return one(G)
    end
    
    if isone(x)
        return x
    else
        return G(n_mod * x.x)
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

Base.one(g::ECGroup{P}) where P <: ECPoint = ECGroup{P}(zero(P))
Base.one(::Type{ECGroup{P}}) where P <: ECPoint = ECGroup{P}(zero(P))

value(g::ECGroup) = value(g.x)

iscompressable(g::ECGroup) = iscompressable(g.x)

struct PGroup{S} <: Group
    g::BigInt

    function PGroup{S}(x_int::Integer; allow_one::Bool=false, skip_validation=false) where S

        x = convert(BigInt, x_int)

        if !allow_one && x == 1
            msg = "Constructing a degenerate element. Use `allow_one` to hide this warning"
            if isstrict()
                error(msg)
            else
                @warn msg
            end
        end
        
        _order = order(PGroup{S})
        _modulus = modulus(PGroup{S})
        if !skip_validation && !isnothing(_order)
            powermod(x, _order, _modulus) == 1 || throw(ArgumentError("Element $x is not an element of prime group with order $_order and modulus $_modulus"))
        end
        
        # Relaxing a little so that group products could happen properly
        0 < x < S.p || throw(ArgumentError("Element $x is not in range of the prime group with modulus $_modulus"))
        new{S}(x)
    end

    #PGroup{S}(x::Integer; skip_validation=false) where S = PGroup{S}(BigInt(x); skip_validation)

    Base.one(::Type{PGroup{S}}) where S = new{S}(1)
end

#Base.one(::Type{PGroup{S}}) where S = PGroup{S}(one)
Base.one(::PGroup{S}) where S = one(PGroup{S})

octet(x::PGroup) = int2octet(value(x), bitlength(modulus(x)))

Base.convert(group::Type{<:PGroup}, element::Vector{UInt8}) = convert(group, octet2int(element))
PGroup{S}(element::Vector{UInt8}) where S = convert(PGroup{S}, element)

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
        end
        return one(G)
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

