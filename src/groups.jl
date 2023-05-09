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

Base.convert(::Type{Vector{G}}, x::Vector) where G <: Group = G[ G <| i for i in x]

Base.rand(prg::PRG, ::Type{G}, N::Integer; nr::Integer = 0) where G <: Group = Vector{G} <| rand(prg, spec(G), N; nr)

Base.ones(x::Vector{G}) where G <: Group = [one(i) for i in x] # need to extend this

struct ECGroup{P<:ECPoint} <: Group
    x::P
end

Base.convert(::Type{ECGroup{P}}, x) where P <: ECPoint = ECGroup{P}(P <| x)
Base.convert(::Type{G}, x::G) where G <: ECGroup = x

Base.:*(x::G, y::G) where G <: ECGroup = G(x.x + y.x)

function Base.:^(x::G, n::Integer) where G <: ECGroup
    n_mod = mod(n, order(G))

    @assert n_mod != 0 "A bad exponent"

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


function Base.show(io::IO, g::G) where G <: ECGroup
    show(io, G)
    print(io, " <| (")
    show(io, gx(g))
    print(io, ", ")
    show(io, gy(g))
    print(io, ")")
end


struct PGroup{S} <: Group
    g::BigInt

    PGroup(p::Integer, q::Integer; name=nothing) = PGroup{static(; p, q, name)}

    function PGroup{S}(x::BigInt) where S
        #@assert 1 < x < S.p "Not in range"
        # Relaxing a little so that group products could happen propeerly
        # But is a presence of 1 in multiplicative products a vulnerability?
        @assert 0 < x < S.p "Not in range" 
        new{S}(x)
    end

    PGroup{S}(x::Integer) where S = PGroup{S}(BigInt(x))
    PGroup{S}(::typeof(one)) where S = new{S}(1)
end


Base.one(::Type{PGroup{S}}) where S = PGroup{S}(one)
Base.one(::PGroup{S}) where S = PGroup{S}(one)

#Base.ones(::         

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

Base.convert(::Type{P}, x::Integer) where P <: PGroup = P(BigInt(x))

Base.isvalid(g::G) where G <: PGroup = value(g) != 1 && powermod(value(g), order(G), modulus(G)) == 1


import Base.*
*(x::G, y::G) where G <: PGroup = G(mod(value(x) * value(y), modulus(G)))


import Base.^
function ^(x::G, n::Integer) where G <: PGroup 

    n_mod = mod(n, order(G))
    @assert n_mod != 0 "A bad exponent"
    #@assert value(x) != 1 "A value 1 is not part of prime group"

    return G(powermod(value(x), n_mod, modulus(G)))
end


import Base.inv
inv(x::G) where G <: PGroup = G(modinv(value(x), modulus(G)))


Base.:(==)(x::G, y::G) where G <: PGroup = x.g == y.g


Base.isless(x::G, y::G) where G <: PGroup = value(x) < value(y)

#Base.convert(::Type{G}, x::Integer) where G <: PGroup = G(x)


# function Base.prod(x::Vector{G}) where G <: PGroup

#     p = modulus(G)1
#     s = value(x[1])
    
#     for i in x[2:end]
#         s = mod(s * value(i), p)
#     end
    
#     return G(s)
# end
