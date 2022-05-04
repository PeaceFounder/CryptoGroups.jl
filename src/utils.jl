macro bin_str(x)
    a = BitVector((i for i in x) .== '1')
    return a
end


tobig(x) = parse(BigInt, bytes2hex(reverse(x)), base=16)

function int2bytes(x::Integer)
    hex = string(x, base=16)
    if mod(length(hex), 2) != 0
        hex = string("0", hex)
    end
    
    return reverse(hex2bytes(hex))
end


struct StaticBigInt{N} <: Integer
    x::NTuple{N, UInt8}
end

# function StaticBigInt(x::Integer)
#     bytes = int2bytes(x)
#     N = length(bytes)
#     sv = NTuple{N, UInt8}(bytes)
#     return StaticBigInt(sv)
# end


function StaticBigInt(x::Integer; n = bitlength(x))
    bytes = int2bytes(x)
    
    N = div(n, 8, RoundUp)

    append!(bytes, UInt8[0 for i in 1:N-length(bytes)])

    #N = length(bytes)

    sv = NTuple{N, UInt8}(bytes)
    return StaticBigInt(sv)
end



Base.convert(::Type{BigInt}, x::StaticBigInt) = tobig(x.x)
Base.BigInt(x::StaticBigInt) = tobig(x.x)


style(x, n) = "\33[1;$(n)m$x\33[0m"

function Base.show(io::IO, a::StaticBigInt)
    show(io, BigInt(a))
end

Base.display(x::StaticBigInt) = show(x)



nbits(x::UInt8) = 8
nbits(x::UInt64) = 64
nbits(x::UInt32) = 32
nbits(x::BigInt) = 64 * x.size

nbits(x::Int64) = 64
nbits(x::Int32) = 32
nbits(x::Int128) = 128

function tobin(a::Integer) 
    s = BitVector(a << -i % Bool for i in 0:(nbits(a)-1))
    n = findlast(x->x==1, s)
    if n == nothing
        n = 0
    end
    return reverse(s[1:n])
end


function tobin(a::Integer, n::Integer) 
    s = BitVector(a << -i % Bool for i in 0:(nbits(a)-1))
    
    if length(s) < n
        resize!(s, n)
        return reverse(s)
    else
        return reverse(s[1:n])
    end
end


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



function uint2string(x::Unsigned)
    bytes = reinterpret(UInt8, [x])

    N = findlast(x->x!=0, bytes)
    
    if N == nothing
        return String("") 
    else
        return String(bytes[1:N]) 
    end
    #trimmed = bytes[1:findlast(x->x!=0, bytes)]
end

function string2uint(x::String)
    bytes = Vector{UInt8}(x)
    padded_bytes = append!(bytes, UInt8[0 for i in 1:16-length(bytes)])
    uint = reinterpret(UInt128, padded_bytes)
    return uint[1]
end
