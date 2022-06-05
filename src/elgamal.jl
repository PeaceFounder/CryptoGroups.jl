struct Enc{T<:Group} 
    pk::T
    g::T
end


(enc::Enc{T})(m::T, r::Integer) where T <: Group = (enc.g^r, m*enc.pk^r) ### Message first?
(enc::Enc)(r::Integer) = (enc.g^r, enc.pk^r)  


a(x::Tuple{T, T}) where T <: Group = x[1]
b(x::Tuple{T, T}) where T <: Group = x[2]


*(x::Tuple{G, G}, y::Tuple{G, G}) where G <: Group = (a(x)*a(y), b(x)*b(y))

(enc::Enc)(e::Tuple{G, G}, r::Integer) where G <: Group = e * enc(r)


struct ElGamal{G <: Group} <: AbstractVector{G}
    a::Vector{G}
    b::Vector{G}

    function ElGamal{G}(a::Vector{G}, b::Vector{G}) where {G <: Group} 
        @assert length(a) == length(b)
        return new(a, b)
    end
end


ElGamal(a::Vector{G}, b::Vector{G}) where G <: Group = ElGamal{G}(a, b)

ElGamal(e::Vector{Tuple{G, G}}) where G <: Group = ElGamal([a(i) for i in e], [b(i) for i in e])

function ElGamal{G}(a::Vector{T}, b::Vector{T}) where {T, G<:Group}
    a′ = convert(Vector{G}, a)
    b′ = convert(Vector{G}, b)

    return ElGamal{G}(a′, b′)
end


a(e::ElGamal) = e.a
b(e::ElGamal) = e.b

Base.getindex(e::ElGamal, i::Int) = (a(e)[i], b(e)[i])
Base.getindex(e::ElGamal{G}, ivec::Vector) where G <: Group = ElGamal{G}(a(e)[ivec], b(e)[ivec])
Base.length(e::ElGamal) = length(a(e))
Base.size(e::ElGamal) = size(a(e))

Base.copy(e::ElGamal{G}) where G <: Group = ElGamal{G}(copy(a(e)), copy(b(e)))
Base.sort(e::ElGamal) = sort!(copy(e))

### NEW METHOD
function Base.setindex!(e::ElGamal{G}, val::Tuple{G, G}, i) where G <: Group
    a(e)[i] = a(val)
    b(e)[i] = b(val)
end


function *(x::ElGamal{G}, y::ElGamal{G}) where G <: Group

    @assert length(x) == length(y)

    a′ = a(x) .* a(y)
    b′ = b(x) .* b(y)

    return ElGamal(a′, b′)
end

function *(x::ElGamal{G}, y::Tuple{G, G}) where G <: Group
    
    a′ = a(x) .* a(y)
    b′ = b(x) .* b(y)

    return ElGamal(a′, b′)
end

*(x::Tuple{G, G}, y::ElGamal{G}) where G <: Group = y * x


import Base: ^

^(x::Tuple{G, G}, k::Integer) where G <: Group = (x[1]^k, x[2]^k)
^(x::ElGamal, k::AbstractVector{<:Integer}) = ElGamal(a(x) .^ k, b(x) .^ k)
^(x::ElGamal, k::Integer) = ElGamal(a(x) .^ k, b(x) .^ k)

Base.broadcasted(::typeof(^), x::ElGamal, y::AbstractVector{<:Integer}) = x^y


(enc::Enc)(e::ElGamal, r::Integer) = enc(r) * e 


function (enc::Enc{G})(m::Vector{G}, r::AbstractVector{<:Integer}) where G <: Group

    a′ = enc.g .^ r
    b′ = m .* (enc.pk .^ r)

    return ElGamal(a′, b′)
end


function (enc::Enc{G})(𝐞::ElGamal{G}, 𝐫::AbstractVector{<:Integer}) where G <: Group
    
    ### Need to improve this 
    r = ElGamal(enc.(𝐫))
    𝐞′ = 𝐞 * r
    
    return 𝐞′    
end


struct Dec
    sk::Integer
end

(dec::Dec)(e::Tuple{G, G}) where G = b(e) * a(e)^(-dec.sk) 

(dec::Dec)(e::ElGamal) = [dec(ei) for ei in e]

Base.isless(x::Tuple{G, G}, y::Tuple{G, G}) where G <: Group = x[1] == y[1] ? x[2] < y[2] : x[1] < y[1]


Base.prod(e::ElGamal{G}) where G <: Group = (prod(e.a), prod(e.b))
