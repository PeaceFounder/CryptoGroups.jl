module Utils

import CryptoPRG: bitlength
import Base.GMP

"""
    @check(ex, msg = nothing)

Drop in replacement for `@assert` macro as it can be eliminated in optimizations.
"""
macro check(ex, msg = nothing)
    return quote
        if !$(esc(ex))
            failed_expr = $(QuoteNode(ex))
            error_msg = $(msg === nothing ? string(ex) : msg)
            file = $(string(__source__.file))
            line = $(__source__.line)
            throw(AssertionError("$error_msg"))
        end
        nothing
    end
end

function jacobi(n::BigInt, k::BigInt)

    # Get the last limb (word) of a BigInt
    get_limb(n::BigInt, i::Int) = unsafe_load(n.d, 1 + i)

    iseven(k) && throw("Argument k=$k should be odd.")
    n = mod(n, k)
    t = 1
    k_1 = true  # k is odd by requirement
    last_limb = get_limb(k, 0)
    k_2 = (last_limb & 2 == 2)  # second bit
    k_3 = (last_limb & 4 == 4)  # third bit
    
    while !iszero(n)
        z = trailing_zeros(n)
        n >>= z
        # After right shift, check bits directly on the last limb
        last_limb = get_limb(n, 0)
        n_2 = (last_limb & 2 == 2)  # second bit
        
        # For k ≡ 3,5 (mod 8) when z is odd
        if isodd(z) && (k_2 ⊻ k_3)
            t = -t
        end
        
        # For quadratic reciprocity when both are ≡ 3 (mod 4)
        if k_2 && n_2
            t = -t
        end
        
        # Save current n's bits for next iteration's k
        k_2 = n_2
        k_3 = (last_limb & 4 == 4)
        n, k = k, n
        n = mod(n, k)
    end
    return isone(k) ? t : 0
end


"""
    int2octet(x::Integer, N::Int = bitlength(x))::Vector{UInt8}

Converts integer `x` into an octet where optional `N` specifies number of allocated bits for the integer encoding. 
"""
function int2octet(x::Integer, N::Int) # We will need a refactor here

    error("This method will be deprecated")
    #@warn "This method will be deprecated"

    k = div(N, 8, RoundUp)
    
    # Pre-allocate the full buffer with zeros
    result = zeros(UInt8, k)
    
    # Convert to hex and ensure even length
    hex = string(x, base=16)
    if mod(length(hex), 2) != 0
        hex = string("0", hex)
    end
    
    # Fill the latter part of result with the converted bytes
    hex2bytes!(view(result, k-div(length(hex),2)+1:k), hex)
    
    return result
end

function int2octet!(buffer::Union{Vector{UInt8}, SubArray{UInt8, 1}}, n::BigInt)

    if iszero(n)
        fill!(buffer, 0)
        return 1
    end
    
    # Calculate number of bytes needed
    nbits = ndigits(n, base=2)
    needed_bytes = cld(nbits + (n < 0 ? 1 : 0), 8)
    
    # Fill leading positions with zeros
    if needed_bytes < length(buffer)
        fill!(@view(buffer[1:length(buffer)-needed_bytes]), 0)
    end
    
    # Export bytes directly to the end of buffer
    _, written = GMP.MPZ.export!(@view(buffer[end-needed_bytes+1:end]), n;
                                 order = 1,    # Big-endian
                                 endian = 0,   # Native endian
                                 nails = 0     # Use all bits
                                 )
    
    return written
end

function int2octet(n::BigInt)

    nbits = ndigits(n, base=2)
    # Add one bit for sign in case of negative numbers
    nbytes = cld(nbits + (n < 0 ? 1 : 0), 8)
    
    # Allocate output buffer
    buffer = Vector{UInt8}(undef, nbytes)
    
    int2octet!(buffer, n)

    return buffer
end

int2octet(n::Integer) = reverse(reinterpret(UInt8, [n]))

function int2octet!(buffer::Union{Vector{UInt8}, SubArray{UInt8, 1}}, n::Integer)

    bytes = int2octet(n)
    copyto!(buffer, @view(bytes[end-length(buffer)-1:end]))

    return
end


"""
    octet2int(x::Vector{UInt8})::BigInt
    octet2int(x::String)::BigInt

Converts a binary octet to a `BigInt`. In case a string is passed it is treated as hexadecimal and is converted with `hex2bytes`. Equivalent in represeentation as `parse(BigInt, bytes2hex(x), base=16)`.
"""
function octet2int(bytes::AbstractVector{UInt8})
    isempty(bytes) && return BigInt(0)
    
    #result = BigInt()
    n_bytes = length(bytes)
    
    # Calculate required limbs
    limb_bytes = sizeof(GMP.Limb)
    nlimbs = cld(n_bytes, limb_bytes)
    
    # Initialize BigInt with pre-allocated size
    result = BigInt(; nbits = nlimbs * GMP.BITS_PER_LIMB)
    
    # Process full limbs first using direct memory access
    if n_bytes >= limb_bytes
        # Create pointer to input bytes for direct memory access
        bytes_ptr = pointer(bytes)
        
        # Process full limbs
        full_limbs = n_bytes ÷ limb_bytes
        for i in 1:full_limbs
            # Calculate offset from end of array
            offset = n_bytes - i * limb_bytes
            # Load limb directly from memory with proper byte order
            limb = unsafe_load(Ptr{GMP.Limb}(bytes_ptr + offset))
            # Handle endianness if needed
            if ENDIAN_BOM == 0x04030201
                limb = ntoh(limb)
            end
            unsafe_store!(result.d, limb, i)
        end
        
        # Handle remaining bytes in the last partial limb
        remaining_bytes = n_bytes % limb_bytes
        if remaining_bytes > 0
            limb = zero(GMP.Limb)
            for j in 1:remaining_bytes
                shift = (remaining_bytes - j) * 8
                limb |= GMP.Limb(bytes[j]) << shift
            end
            unsafe_store!(result.d, limb, nlimbs)
        end
    else
        # Handle case when input is smaller than a limb
        limb = zero(GMP.Limb)
        for j in 1:n_bytes
            shift = (n_bytes - j) * 8
            limb |= GMP.Limb(bytes[j]) << shift
        end
        unsafe_store!(result.d, limb, 1)
    end
    
    # Set the correct size (number of non-zero limbs)
    actual_limbs = nlimbs
    while actual_limbs > 0 && unsafe_load(result.d, actual_limbs) == 0
        actual_limbs -= 1
    end
    result.size = actual_limbs
    
    return result
end

#octet2int(x::Vector{UInt8}) = parse(BigInt, bytes2hex(x), base=16)

octet2int(x::String) = octet2int(hex2bytes(x))
octet2int(x::NTuple{N, UInt8}) where N = octet2int([x...])

"""
    octet2bits(x::Vector{UInt8}[, N::Int])::BitVector

Converts an octet ot a bitvector with an optionally specified length `N`.
"""    
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

"""
    bits2octet(x::BitVector)::Vector{UInt8}

Converts a bitvector to an octet.
"""
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

"""
    @bin_str(x)::BitVector

Converts binary representation to `BitVector`:
```julia
bin"1010" == BitVector([1, 0, 1, 0])
```
"""
macro bin_str(x)
    a = BitVector((i for i in x) .== '1')
    return a
end

"""
    dynamic(x)

Converts static type to it's closest base type. Inverse of `static` method.
"""
dynamic(x) = x 

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

dynamic(x::StaticBigInt) = convert(BigInt, x)

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

"""
    @hex_str(x)::Vector{UInt8}

A convinience macro for converting a hex to byte vector:
```julia
hex"AAEF BBEC" == UInt8[0xaa, 0xef, 0xbb, 0xec]
```
"""
macro hex_str(x)
    return _hex2bytes(x)
end

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

dynamic(x::StaticBitVector) = convert(BitVector, x)

function Base.show(io::IO, x::StaticBitVector) 
    print(io, "static(")
    print(io, "bin\"")
    print(io, join(i ? "1" : "0" for i in convert(BitVector, x)))
    print(io, "\"")
    print(")")
end

Base.length(x::StaticBitVector) = x.len

########################

"""
    static(x)
    static(; kwargs...)

If `x` is not bitstype converts to a bitstype representation. Keyword arguments are used to construct a named tuple with every value
being converted with static. Can be converted back to original type with `dynamic` method.

# Example

```julia
# Single argument case
x = BigInt(23)
isbitstype(typeof(x)) == false
isbitstype(typeof(static(x))) == true
x == dynamic(static(x))

# NamedTuple construction

nt = static(; x = BigInt(23), y = 51, z = :name)
isbittstype(typeof(nt)) == true
nt.x == BigInt(23) # accessor methods do conversion with dynamic
dynamic(nt) # ordinary named tuple with every value made dynamic
```

"""
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
dynamic(x::StaticSymbol) = convert(Symbol, x)

Base.show(io::IO, x::StaticSymbol) = show(io, dynamic(x))

struct StaticNamedTuple{T, V}
    args::NamedTuple{T, V}
end

StaticNamedTuple(; kwargs...) = StaticNamedTuple(NamedTuple((key, static(value)) for (key, value) in pairs(kwargs)))

Base.propertynames(x::StaticNamedTuple) = propertynames(getfield(x, :args))
Base.getproperty(x::StaticNamedTuple, sym::Symbol) = dynamic(getfield(getfield(x, :args), sym))

static(; kwargs...) = StaticNamedTuple(; kwargs...)
dynamic(x::StaticNamedTuple) = NamedTuple((key, dynamic(value)) for (key, value) in pairs(getfield(x, :args)))

function Base.show(io::IO, x::StaticNamedTuple)

    print(io, "static")
    show(io, dynamic(x))

    return
end





export static, dynamic, @bin_str, @hex_str, octet2int, int2octet, int2octet!, octet2bits, bits2octet

end
