using .Curves: ECPoint, gx, gy
using .Specs: PRG


abstract type Group end

Base.broadcasted(f::Function, x::Group, y::AbstractVector{<:Integer}) = f.((x for i in 1:length(y)), y)

Base.broadcasted(::typeof(*), x::G, y::Vector{G}) where G <: Group = (x for i in 1:length(y)) .* y 
Base.broadcasted(::typeof(*), x::Vector{G}, y::G) where G <: Group =  x .* (y for i in 1:length(x))

order(x::G) where G <: Group = order(G)

import Base./
/(x::G, y::G) where G <: Group = x * inv(y)

name(x::G) where G <: Group = name(G)

Base.rand(prg::PRG, ::Type{G}, N::Integer; nr::Integer = 0) where G <: Group = convert(Vector{G}, rand(prg, spec(G), N; nr))

Base.ones(x::Vector{G}) where G <: Group = [one(i) for i in x] # need to extend this

struct ECGroup{P<:ECPoint} <: Group
    x::P
end


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
    end
    
    if isone(x)
        return x
    else
        return G(n_mod * x.x)
    end
end

Base.inv(g::G) where G <: ECGroup = g^(order(G) - 1)


order(::Type{ECGroup{P}}) where P = order(P)

Base.isvalid(x::ECGroup) = isvalid(x.x)

Base.:(==)(x::G, y::G) where G <: ECGroup = x.x == y.x

modulus(::Type{ECGroup{P}}) where P <: ECPoint = modulus(P)

name(::Type{ECGroup{P}}) where P <: ECPoint = name(P)

Base.isless(x::G, y::G) where G <: ECGroup = isless(x.x, y.x)

Curves.gx(g::ECGroup) = gx(g.x)
Curves.gy(g::ECGroup) = gy(g.x)

Base.one(g::ECGroup{P}) where P <: ECPoint = ECGroup{P}(zero(P))
Base.one(::Type{ECGroup{P}}) where P <: ECPoint = ECGroup{P}(zero(P))

# function Base.show(io::IO, g::G) where G <: ECGroup
#     show(io, G)
#     print(io, " <| (")
#     show(io, gx(g))
#     print(io, ", ")
#     show(io, gy(g))
#     print(io, ")")
# end


struct PGroup{S} <: Group
    g::BigInt

    PGroup(p::Integer, q::Integer; name=nothing) = PGroup{static(; p, q, name)}

    function PGroup{S}(x::BigInt; allow_one::Bool=false) where S

        if !allow_one && x == 1
            msg = "Constructing a degenerate element. Use `allow_one` to hide this warning"
            if isstrict()
                error(msg)
            else
                @warn msg
            end
        end

        # Relaxing a little so that group products could happen propeerly
        @assert 0 < x < S.p "Not in range" 
        new{S}(x)
    end

    PGroup{S}(x::Integer) where S = PGroup{S}(BigInt(x))
    PGroup{S}(::typeof(one)) where S = new{S}(1)
end


Base.one(::Type{PGroup{S}}) where S = PGroup{S}(one)
Base.one(::PGroup{S}) where S = PGroup{S}(one)

octet(x::PGroup) = int2octet(value(x), bitlength(modulus(x)))

Base.convert(group::Type{<:PGroup}, element::Vector{UInt8}) = convert(group, octet2int(element))
PGroup{S}(element::Vector{UInt8}) where S = convert(PGroup{S}, element)


modulus(::Type{PGroup{S}}) where S = BigInt(S.p)
modulus(::G) where G <: PGroup = modulus(G)

order(::Type{PGroup{S}}) where S = BigInt(S.q)

name(::Type{PGroup}) = nothing

name(::Type{PGroup{S}}) where S =  isnothing(S.name) ? nothing : convert(Symbol, S.name)
    

Base.show(io::IO, g::PGroup) = print(io, value(g))


Base.show(io::IO, ::Type{PGroup}) = print(io, "PGroup")

#dublicate at Fields. 
function trimnumber(x::String)
    if length(x) < 30
        return x
    else
        return x[1:10] * "..." * x[end-10:end]
    end
end

trimnumber(x::Integer) = trimnumber(string(x))

groupstr(m) = "𝐙/($(trimnumber(m)))"


function Base.show(io::IO, ::Type{G}) where G <: PGroup
    
    if name(G) == nothing
        print(io, groupstr(modulus(G)))
    else
        print(io, name(G))
    end
end


function Base.display(::Type{G}) where G <: PGroup
    show(G)
    print(" (order = $(order(G)))")
end


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
        return  G(mod(value(x) * value(y), modulus(G)))
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
    end

    if isone(x)
        return x
    else
        return G(powermod(value(x), n_mod, modulus(G)))
    end
end


import Base.inv
inv(x::G) where G <: PGroup = G(modinv(value(x), modulus(G)))


Base.:(==)(x::G, y::G) where G <: PGroup = x.g == y.g


Base.isless(x::G, y::G) where G <: PGroup = value(x) < value(y)

