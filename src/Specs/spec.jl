abstract type Spec end

import Base: @kwdef

#_parse_seed(seed::String) = hex2bytes(seed)
#_parse_seed(seed::Vector{UInt8}) = seed

# Dublicate present in utlis.jl
function _hex2bytes(x::String)
    
    normalized = join(split(x, " "), "") 

    N = length(normalized)
    
    if mod(N, 2) != 0
        normalized = "0" * normalized
    end
    
    return hex2bytes(normalized)
end

function bytes2bits(x::Vector{UInt8})
    bv = BitVector(u << -i % Bool for u in x for i in 7:-1:0)
    return bv
end

hex2bits(x::String) = bytes2bits(_hex2bytes(x))



_parse_bits(x::String, N::Int) = hex2bits(x)[end - N + 1:end]
_parse_bits(x::BitVector, m::Int) = x
_parse_bits(x::Nothing, m::Int) = x


_parse_int(x::String) = parse(BigInt, x, base=16)
_parse_int(x::Integer) = BigInt(x)
_parse_int(::Nothing) = nothing


@kwdef struct ECP <: Spec
    p::BigInt
    n::Union{BigInt, Nothing} = nothing
    a::BigInt = -3
    b::BigInt
    Gx::Union{BigInt, Nothing} = nothing
    Gy::Union{BigInt, Nothing} = nothing

    function ECP(p, n, a, b, Gx, Gy)
        
        _a = mod(_parse_int(a), p) # taking mod as conventually a=-3
        _b = _parse_int(b)
        _Gx = _parse_int(Gx)
        _Gy = _parse_int(Gy)

        #return ECP(p, n, _a, _b, _Gx, _Gy)
        return new(p, n, _a, _b, _Gx, _Gy)
    end
end

order(curve::ECP) = curve.n
generator(curve::ECP) = (curve.Gx, curve.Gy)

modulus(curve::ECP) = curve.p


# I could always add a field for equation to be used

Base.:(==)(x::ECP, y::ECP) = x.p == y.p && x.n == y.n && x.a == y.a && x.b == y.b && x.Gx == y.Gx && x.Gy == y.Gy



abstract type BinaryBasis end

struct PB <: BinaryBasis
    f::Vector{Int}
end

bitlength(x::PB) = maximum(x.f) #- 1

struct GNB <: BinaryBasis
    m::Int
    T::Int
end

bitlength(x::GNB) = x.m


struct EC2N{B<:BinaryBasis} <: Spec
    basis::B
    n::BigInt
    a::BitVector
    b::BitVector
    Gx::Union{BitVector, Nothing}
    Gy::Union{BitVector, Nothing}

    function EC2N(basis::B, n::BigInt, a::BitVector, b::BitVector, Gx::BitVector, Gy::BitVector) where B <: BinaryBasis
        @assert bitlength(basis) == length(a) == length(b) == length(Gx) == length(Gy) 

        new{B}(basis, n, a, b, Gx, Gy)
    end
    
    function EC2N(basis::B, n::BigInt, a::BitVector, b::BitVector) where B <: BinaryBasis
        @assert bitlength(basis) == length(a) == length(b) 
        new{B}(basis, n, a, b, nothing, nothing)
    end
end


_parse_bits(x::BitVector, basis::BinaryBasis) = x
_parse_bits(x::String, basis::BinaryBasis) = _parse_bits(x, bitlength(basis))


function _int2bits_gnb(x::Int, m::Int)
    if x == 0
        return BitVector(0 for i in 1:m)
    elseif x == 1
        return BitVector(1 for i in 1:m)
    else
        error("Conversion of $x not possible")
    end
end

_parse_bits(x::Int, basis::GNB) = _int2bits_gnb(x, bitlength(basis))



function _int2bits_pb(x::Int, m::Int)
    if x == 0
        return BitVector(0 for i in 1:m)
    elseif x == 1
        return BitVector(((0 for i in 1:m-1)..., 1))
    else
        error("Conversion of $x not possible")
    end
end

_parse_bits(x::Int, basis::PB) = _int2bits_pb(x, bitlength(basis))



# Basis could be nonoptional argument here
function EC2N(basis::BinaryBasis; n, a, b, Gx=nothing, Gy=nothing) 

    _n = convert(BigInt, n)
    _a = _parse_bits(a, basis)
    _b = _parse_bits(b, basis)
    _Gx = _parse_bits(Gx, bitlength(basis))
    _Gy = _parse_bits(Gy, bitlength(basis))

    if isnothing(_Gx) && isnothing(_Gy)
        return EC2N(basis, _n, _a, _b)
    elseif !isnothing(_Gx) && !isnothing(_Gy)
        return EC2N(basis, _n, _a, _b, _Gx, _Gy)
    else
        error("Incompatible input.")
    end
end


function EC2N{PB}(; f, kwargs...)
    basis = PB(f)
    return EC2N(basis; kwargs...)
end

function EC2N{GNB}(; m, T, kwargs...)
    basis = GNB(m, T)
    return EC2N(basis; kwargs...)
end

function EC2N(; kwargs...)

    if :f in kwargs
        return EC2N{PB}(kwargs...)
    elseif :m in kwargs && :T in kwargs
        return EC2N{GNB}(kwargs...)
    else
        error("Dispatch can not be infered from passed keyword arguments")
    end
end


order(curve::EC2N) = curve.n
generator(curve::EC2N) = (curve.Gx, curve.Gy)

a(curve::EC2N) = curve.a
b(curve::EC2N) = curve.b


struct Koblitz{B<:BinaryBasis} <: Spec
    bec::EC2N{B}

    function Koblitz{B}(; a, kwargs...) where B <: BinaryBasis

        b = 1
        @assert a in [0, 1]

        bec =  EC2N{B}(; a, b, kwargs...)

        return new{B}(bec)
    end

    function Koblitz(; a, kwargs...)

        b = 1
        @assert a in [0, 1]

        bec =  EC2N(; a, b, kwargs...)

        B = typeof(bec.basis)

        return new{B}(bec)
    end
end


order(curve::Koblitz) = order(curve.bec)
generator(curve::Koblitz) = generator(curve.bec)

a(curve::Koblitz) = a(curve.bec)
b(curve::Koblitz) = b(curve.bec)


### MODP spec 

tobint(x::String) = parse(BigInt, x, base=16)
tobint(x::BigInt) = x
tobint(x::Integer) = BigInt(x)
tobint(x::Nothing) = nothing

struct MODP <: Spec
    p::BigInt
    g::BigInt
    q::Union{BigInt, Nothing}
    
    MODP(p, g, q) = new(tobint(p), tobint(g), tobint(q))    
    MODP(p, g) = MODP(p, g, nothing)
end


generator(spec::MODP) = spec.g




# abstract type EC2N <: Spec end

# order(curve::EC2N) = curve.n
# generator(curve::EC2N) = (curve.Gx, curve.Gy)

# a(curve::EC2N) = curve.a
# b(curve::EC2N) = curve.b 


# @kwdef struct Koblitz_GNB <: EC2N
#     m::Int
#     T::Int
#     n::BigInt
#     a::Int
#     Gx::BitVector
#     Gy::BitVector

#     function Koblitz_GNB(m::Int, T::Int, n::BigInt, a::Int, Gx::BitVector, Gy::BitVector)
        
#         @assert a in [0, 1]
#         @assert length(Gx) == length(Gy) == m
        
#         new(m, T, n, a, Gx, Gy)
#     end
# end

# function Koblitz_GNB(m, T, n, a, Gx, Gy)

#     _m = convert(Int, m)
#     _T = convert(Int, T)
#     _a = convert(Int, a)
#     _n = convert(BigInt, n)
#     _Gx = _parse_bits(Gx, _m)
#     _Gy = _parse_bits(Gy, _m)

#     return Koblitz_GNB(_m, _T, _n, _a, _Gx, _Gy)
# end


# function b(curve::Koblitz_GNB)
    
#     (; m) = curve

#     _b = 1
#     b′ = _parse_bits_gnb(_b, m)

#     return b′
# end


# @kwdef struct Koblitz_PB <: EC2N
#     f::Vector{Int}
#     n::BigInt
#     a::Int
#     Gx::BitVector
#     Gy::BitVector

#     function Koblitz_PB(f::Vector{Int}, n::BigInt, a::Int, Gx::BitVector, Gy::BitVector)
        
#         m = maximum(f)

#         @assert a in [0, 1]
#         @assert length(Gx) == length(Gy) == m

#         new(f, n, a, Gx, Gy)
#     end
# end

# function Koblitz_PB(f, n, a, Gx, Gy)

#     _f = convert(Vector{Int}, f)

#     m = maximum(_f)

#     _n = convert(BigInt, n)
#     _a = convert(Int, a)
#     _Gx = _parse_bits(Gx, m)
#     _Gy = _parse_bits(Gy, m)

#     return Koblitz_PB(_f, _n, _a, _Gx, _Gy)
# end


# function b(curve::Koblitz_PB) 

#     (; f) = curve

#     m = maximum(f)

#     _b = 1
#     b′ = _parse_bits_pb(_b, m)
    
#     return b′
# end




# #abstract type BinaryCurveSpec <: EllipticCurveSpec end

# @kwdef struct BEC_GNB <: EC2N  #BinaryCurveSpec
#     m::Int
#     T::Int
#     n::BigInt
#     a::BitVector
#     b::BitVector
#     Gx::BitVector
#     Gy::BitVector

#     function BEC_GNB(m::Int, T::Int, n::BigInt, a::BitVector, b::BitVector, Gx::BitVector, Gy::BitVector)
        
#         @assert m == length(a) == length(b) == length(Gx) == length(Gy)

#         new(m, T, n, a, b, Gx, Gy)
#     end
# end

# function BEC_GNB(m, T, n, a, b, Gx, Gy)
    
#     _m = convert(Int, m)
#     _T = convert(Int, T)
#     _n = convert(BigInt, n)
#     _a = _parse_bits_gnb(a, _m)
#     _b = _parse_bits_gnb(b, _m)
#     _Gx = _parse_bits(Gx, _m)
#     _Gy = _parse_bits(Gy, _m)

#     return BEC_GNB(_m, _T, _n, _a, _b, _Gx, _Gy)
# end


# function BEC_GNB(curve::Koblitz_GNB)
    
#     (; m, T, n, a, Gx, Gy) = curve

#     _b = b(curve)

#     return BEC_GNB(m, T, n, a, _b, Gx, Gy)
# end





# @kwdef struct BEC_PB <: EC2N  #BinaryCurveSpec
#     f::Vector{Int}
#     n::BigInt
#     a::BitVector   
#     b::BitVector
#     Gx::BitVector
#     Gy::BitVector

#     function BEC_PB(f::Vector{Int}, n::BigInt, a::BitVector, b::BitVector, Gx::BitVector, Gy::BitVector)
        
#         m = maximum(f)

#         @assert length(a) == length(b) == length(Gx) == length(Gy) == m

#         new(f, n, a, b, Gx, Gy)
#     end
# end


# function BEC_PB(f, n, a, b, Gx, Gy)
    
#     _f = convert(Vector{Int}, f)
#     _m = maximum(_f)
#     _n = convert(BigInt, n)
    
#     _a = _parse_bits_pb(a, _m)
#     _b = _parse_bits_pb(b, _m)
#     _Gx = _parse_bits(Gx, _m)
#     _Gy = _parse_bits(Gy, _m)

#     return BEC_PB(_f, _n, _a, _b, _Gx, _Gy)
# end



# function BEC_PB(curve::Koblitz_PB)
    
#     (; f, n, a, Gx, Gy) = curve

#     _b = b(curve)

#     return BEC_PB(f, n, a, _b, Gx, Gy)
# end


