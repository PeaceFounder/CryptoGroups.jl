module Utils

import CryptoPRG: bitlength

function int2octet(x::Integer)

    hex = string(x, base=16)
    if mod(length(hex), 2) != 0
        hex = string("0", hex)
    end
    
    return hex2bytes(hex)
end


function int2octet(x::Integer, N::Int)
    k = div(N, 8, RoundUp)


    bytes = int2octet(x)

    pad = UInt8[0 for i in 1:(k - length(bytes))]

    return UInt8[pad..., bytes...]
end


octet2int(x::Vector{UInt8}) = parse(BigInt, bytes2hex(x), base=16)
octet2int(x::String) = octet2int(hex2bytes(x))

octet2int(x::NTuple{N, UInt8}) where N = octet2int([x...])
    
function octet2bits(x::Vector{UInt8})
    bv = BitVector(u << -i % Bool for u in x for i in 7:-1:0)
    return bv
end

octet2bits(x::Vector{UInt8}, N::Int) = octet2bits(x)[end - N + 1:end]
octet2bits(x::String, N::Int) = octet2bits(hex2bytes(x), N)


function bits2uint8(x::BitVector)

    s = UInt8(0)

    for (i, j) in enumerate(reverse(x))

        if j == 1
            s += UInt8(2)^(i - 1)
        end
    end

    return s
end

function bits2octet(_x::BitVector)
    
    x = copy(_x)

    if mod(length(x), 8) != 0

        padding = BitVector(0 for i in 1:(8 - mod(length(x), 8)))
        prepend!(x, padding)

    end

    # For now assuming that x is in the length of octets

    N = div(length(x), 8, RoundUp)

    b = reshape(x, 8, N)

    bytes = UInt8[bits2uint8(b[:, i]) for i in 1:N]

    return bytes
end


macro bin_str(x)
    a = BitVector((i for i in x) .== '1')
    return a
end

struct StaticBigInt{N} <: Integer
    x::NTuple{N, UInt8}
end

function StaticBigInt(x::Integer; n = bitlength(x))
    bytes = int2octet(x)
    
    N = div(n, 8, RoundUp)

    append!(bytes, UInt8[0 for i in 1:N-length(bytes)])

    sv = NTuple{N, UInt8}(bytes)
    return StaticBigInt(sv)
end


Base.convert(::Type{BigInt}, x::StaticBigInt) = octet2int(x.x)
Base.convert(::Type{Integer}, x::StaticBigInt) = convert(BigInt, x)

Base.BigInt(x::StaticBigInt) = octet2int(x.x)

style(x, n) = "\33[1;$(n)m$x\33[0m"

function Base.show(io::IO, a::StaticBigInt)
    show(io, BigInt(a))
end

Base.display(x::StaticBigInt) = show(x)


#######################

function _hex2bytes(x::String)
    
    normalized = join(split(x, " "), "") 

    N = length(normalized)
    
    if mod(N, 2) != 0
        normalized = "0" * normalized
    end
    
    return hex2bytes(normalized)
end

macro hex_str(x)
    return _hex2bytes(x)
end


# function bytes2bits(x::Vector{UInt8})
#     bv = BitVector(u << -i % Bool for u in x for i in 7:-1:0)
#     return bv
# end

hex2bits(x::String) = octet2bits(_hex2bytes(x))

tobits(x::String) = hex2bits(x)
tobits(x::Vector{UInt8}) = bytes2bits(x)
tobits(x) = convert(BitVector, x)


nbits(x::Integer) = sizeof(x) * 8
nbits(x::BigInt) = x.size * 64

function tobits(a::Integer) # I could still associate keyword arguments
    s = BitVector(a << -i % Bool for i in 0:(nbits(a)-1))
    n = findlast(x->x==1, s)
    if n == nothing
        n = 0
    end
    return reverse(s[1:n])
end


function tobits(a::Integer, n::Integer) 
    s = BitVector(a << -i % Bool for i in 0:(nbits(a)-1))
    
    if length(s) < n
        resize!(s, n)
        return reverse(s)
    else
        return reverse(s[1:n])
    end
end


#####################


struct StaticBitVector{N}
    len::Int
    chunks::NTuple{N, UInt64}
end

StaticBitVector(x::BitVector) = StaticBitVector(x.len, Tuple(x.chunks))

function Base.convert(::Type{BitVector}, x::StaticBitVector)

    a = BitVector(undef, x.len)
    a.chunks = UInt64[x.chunks...]
    
    return a
end


Base.show(io::IO, x::StaticBitVector) = print(io, join(i ? "1" : "0" for i in convert(BitVector, x)))


modinv(s, q) = mod(gcdx(s, q)[2], q)

Base.length(x::StaticBitVector) = x.len

########################

static(x::BitVector) = StaticBitVector(x)
static(x::BigInt) = StaticBigInt(x)
static(x::Integer) = x # BigInt is the only exception
static(::Nothing) = nothing


struct StaticSymbol
    x::UInt128

    function StaticSymbol(s::Symbol)
        x = string(s)
        bytes = Vector{UInt8}(x)
        padded_bytes = append!(bytes, UInt8[0 for i in 1:16-length(bytes)])
        uint = reinterpret(UInt128, padded_bytes)
        return new(uint[1])
    end
end

function Base.convert(::Type{Symbol}, s::StaticSymbol) 
    
    x = s.x

    bytes = reinterpret(UInt8, [x])

    N = findlast(x->x!=0, bytes)
    
    if N == nothing
        error("Nothing encoded")
    else
        return Symbol(String(bytes[1:N]))
    end
end

static(x::Symbol) = StaticSymbol(x)

# Goeas through the keyword arguments and returns a coresponding namedtuple where on each argument static is being called
static(; kwargs...) = NamedTuple((key, static(value)) for (key, value) in pairs(kwargs))


export static, @bin_str, @hex_str, octet2int, int2octet, octet2bits, bits2octet

end
