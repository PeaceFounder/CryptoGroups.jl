abstract type Group end

Base.broadcasted(f::Function, x::Group, y::AbstractVector{<:Integer}) = f.((x for i in 1:length(y)), y)

Base.broadcasted(::typeof(*), x::G, y::Vector{G}) where G <: Group = (x for i in 1:length(y)) .* y 
Base.broadcasted(::typeof(*), x::Vector{G}, y::G) where G <: Group =  x .* (y for i in 1:length(x))


struct ECGroup{P<:ECPoint} <: Group
    x::P
end

<|(::Type{ECGroup{P}}, x) where P <: ECPoint = ECGroup{P}(P <| x)

Base.:*(x::G, y::G) where G <: ECGroup = G(x.x + y.x)

function Base.:^(x::G, n::Integer) where G <: ECGroup
    if n > 0
        return G(n * x.x)
    elseif n < 0
        m = order(G) + n #+ 1
        return G(m * x.x)
    elseif mod(n, order(G)) == 0
        error("Power leads to a point at infinity")
    end
end

Base.inv(g::G) where G <: ECGroup = g^(order(G) - 1)


order(::Type{ECGroup{P}}) where P = order(P)
order(x::G) where G <: ECGroup = order(G)

validate(x::ECGroup) = validate(x.x)

Base.:(==)(x::G, y::G) where G <: ECGroup = x.x == y.x


gx(p::ECGroup) = gx(p.x)
gy(p::ECGroup) = gy(p.x)

Base.isless(x::G, y::G) where G <: ECGroup = gx(x) == gx(y) ? gx(x) < gx(y) : gy(x) < gy(y)


# Static fields for PGroup
struct static_PGroup{N} ### Type parameter is essential to ensure it to be bitstype
    p::StaticBigInt{N}
    q::StaticBigInt{N} # for the cases where order is unkown I could have a 0 here. 
    name::UInt128

    static_PGroup(p::StaticBigInt{N}, q::StaticBigInt{N}, name::Symbol) where N = new{N}(p, q, string2uint(string(name)))
    function static_PGroup(p::Integer, q::Integer, name::Symbol) 
        # Need to resolve padding to make them equal

        n = bitlength(p)

        _p = StaticBigInt(p; n)
        _q = StaticBigInt(q; n)
        #new(_p, _q, name)
        static_PGroup(_p, _q, name)
    end

    static_PGroup(p::Integer, q::Integer) where N = static_PGroup(p, q, Symbol(""))
end

struct PGroup{S} <: Group
    g::BigInt

    function PGroup{S}(x::BigInt) where S
        @assert 1 < x < S.p "Not in range"
        new{S}(x)
    end

    PGroup{S}(x::Integer) where S = PGroup{S}(BigInt(x))
end


specialize(::Type{PGroup}, p::Integer, q::Integer) = PGroup{static_PGroup(p, q)}
specialize(::Type{PGroup}, q::Integer) = specialize(PGroup, 2*q + 1, q) # Safe prime constructor

specialize(::Type{PGroup}, p::Integer, q::Integer, name::Symbol) = PGroup{static_PGroup(p, q, name)}
specialize(::Type{PGroup}, q::Integer, name::Symbol) = specialize(PGroup, 2*q + 1, q, name) # Safe prime constructor


modulus(::Type{PGroup{S}}) where S = BigInt(S.p)
order(::Type{PGroup{S}}) where S = BigInt(S.q)


name(::Type{PGroup}) = nothing

function name(::Type{PGroup{S}}) where S
    
    str = uint2string(S.name)

    if str == ""
        return nothing
    else
        return Symbol(str)
    end
end


Base.show(io::IO, g::PGroup) = print(io, value(g))


Base.show(io::IO, ::Type{PGroup}) = print(io, "PGroup")

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

<|(::Type{P}, x::Integer) where P <: PGroup = P(BigInt(x))


validate(g::G) where G <: PGroup = value(g) != 1 && powermod(value(g), order(G), modulus(G)) == 1


import Base.*
*(x::G, y::G) where G <: PGroup = G(mod(value(x) * value(y), modulus(G)))


import Base.^
^(x::G, n::Integer) where G <: PGroup = G(powermod(value(x), n, modulus(G)))


import Base.inv
inv(x::G) where G <: PGroup = G(modinv(value(x), modulus(G)))

import Base./
/(x::G, y::G) where G <: PGroup = x * inv(y)


Base.:(==)(x::G, y::G) where G <: PGroup = x.g == y.g


Base.isless(x::G, y::G) where G <: PGroup = value(x) < value(y)

Base.convert(::Type{G}, x::Integer) where G <: PGroup = G(x)
