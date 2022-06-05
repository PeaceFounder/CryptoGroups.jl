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
Base.convert(::Type{Integer}, x::StaticBigInt) = convert(BigInt, x)


Base.BigInt(x::StaticBigInt) = tobig(x.x)


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

function bytes2bits(x::Vector{UInt8})
    bv = BitVector(u << -i % Bool for u in x for i in 7:-1:0)
    return bv
end


hex2bits(x::String) = bytes2bits(_hex2bytes(x))


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
        #return String("") 
        error("Nothing encoded")
        #return nothing
    else
        #return String(bytes[1:N]) 
        return Symbol(String(bytes[1:N]))
    end
    #trimmed = bytes[1:findlast(x->x!=0, bytes)]
end

static(x::Symbol) = StaticSymbol(x)


# Goeas through the keyword arguments and returns a coresponding namedtuple where on each argument static is being called
static(; kwargs...) = NamedTuple((key, static(value)) for (key, value) in pairs(kwargs))

