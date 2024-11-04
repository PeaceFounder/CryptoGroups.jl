# Implements point compression. Currently it is tied to AffinePoint but should be extended for generic case
using ..Fields: modulus, value, bitlength, int2octet!, octet2int, octet2bits, bits2octet, tobits
using CryptoUtils: sqrt_mod_prime
import ..Fields: octet

function decompress_weierstrass(x::BigInt, ỹ::Bool, (a, b)::Tuple{BigInt, BigInt}, p::BigInt)

    α = mod(x^3 + a * x + b, p)

    β = sqrt_mod_prime(α, p)

    # Rightmost bit of β. Perhaps this works:
    # Alternativelly I could do circshift
    if mod(β, 2) == ỹ
        y = β
    else
        y = p - β
    end

    return y
end

function decompress_weierstrass(x::BigInt, ỹ::Bool, (a, b)::Tuple{F, F}) where F <: PrimeField

    p = modulus(F)
    
    return decompress_weierstrass(x, ỹ, (value(a), value(b)), p)
end

function Base.convert(::Type{P}, po::AbstractVector{UInt8}) where P <: AbstractPoint

    pc = po[1]

    if pc == 0
        return zero(P)
    elseif pc == 2 
        ỹ = false
        x = po[2:end]
        return P(x, ỹ)
    elseif pc == 3
        ỹ = true
        x = po[2:end]
        return P(x, ỹ)
    elseif pc == 4 || pc == 6 || pc == 7
        
        l = div(length(po) - 1, 2)
        x = po[2:(l+1)]
        y = po[(l+2):(2l+1)]

        if pc == 4
            return P(x, y)
        elseif pc == 6
            return P(x, y, false)
        elseif pc == 7
            return P(x, y, true)
        end
    else
        error("Wrong PC: $pc")
    end
    
end

function (::Type{P})(x_octet::Vector{UInt8}, y_octet::Vector{UInt8}, ỹ::Bool) where P <: AffinePoint{<:BinaryCurve}
    
    p = P(x_octet, y_octet)

    if isstrict()
        @warn "Sign bit is ignored."
    end

    return p
end

function (::Type{P})(x_octet::Vector{UInt8}, y_octet::Vector{UInt8}, ỹ::Bool) where P <: AffinePoint{<:Weierstrass}

    p = P(x_octet, y_octet)

    (x, y) = convert(Tuple{BigInt, BigInt}, p) #value(p) # Perhaps convert method could be used instead?

    y′ = decompress_weierstrass(x, ỹ, (a(P), b(P)))

    if !(y == y′)
        @warn "Contradictory compression."
    end

    return p
end


function (::Type{P})(x_octet::Vector{UInt8}, ỹ::Bool) where P <: AffinePoint{<:Weierstrass}

    x = octet2int(x_octet)
    y = decompress_weierstrass(x, ỹ, (a(P), b(P)))

    return P(x, y)
end

(::Type{ECPoint{P, S}})(po::Vector{UInt8}) where {P <: AffinePoint, S} = ECPoint{P, S}(P(po))


function _octet(x::BigInt, y::BigInt, N::Int; mode::Symbol = :uncompressed)
    if iszero(x) && iszero(y)
        return UInt8[0]
    end
    
    nbytes = cld(N, 8)
    # Preallocate single output buffer
    # Size is 1 (mode byte) + nbytes (x) + nbytes (y) for the largest case
    out = Vector{UInt8}(undef, 2*nbytes + 1)
    
    if mode == :uncompressed
        out[1] = 0x04
        # Write x and y directly into the output buffer at appropriate offsets
        int2octet!(@view(out[2:nbytes+1]), x)
        int2octet!(@view(out[nbytes+2:2nbytes+1]), y)
        written = 2*nbytes + 1
    elseif mode in [:compressed, :hybrid]
        ỹ = mod(y, 2) % Bool
        
        if mode == :compressed
            out[1] = ỹ ? 0x03 : 0x02
            int2octet!(@view(out[2:nbytes+1]), x)
            written = nbytes + 1
        else # mode == :hybrid
            out[1] = ỹ ? 0x07 : 0x06
            int2octet!(@view(out[2:nbytes+1]), x)
            int2octet!(@view(out[nbytes+2:2nbytes+1]), y)
            written = 2*nbytes + 1
        end
    else
        error("Unrecognized mode $mode")
    end
    
    # Return only the portion of the buffer that was written to
    return resize!(out, written)
end


_octet(x::F, y::F; mode::Symbol = :uncompressed) where F <: PrimeField = _octet(value(x), value(y), bitlength(modulus(F)); mode)

function _octet(x::F, y::F; mode::Symbol = :uncompressed) where F <: BinaryField
    if iszero(x) && iszero(y)
        return UInt8[0]
    end
    
    # Get octets for x and y
    _x = octet(x)
    _y = octet(y)
    
    # Calculate result size and allocate buffer
    nbytes_x = length(_x)
    nbytes_y = length(_y)
    
    # Determine output size based on mode
    outsize = if mode == :compressed
        nbytes_x + 1  # 1 byte for header + x
    else # :uncompressed or :hybrid
        nbytes_x + nbytes_y + 1  # 1 byte for header + x + y
    end
    
    out = Vector{UInt8}(undef, outsize)
    
    if mode == :uncompressed
        # [4; x; y]
        out[1] = 0x04
        copyto!(out, 2, _x, 1, nbytes_x)
        copyto!(out, nbytes_x + 2, _y, 1, nbytes_y)
        
    elseif mode in [:compressed, :hybrid]
        if isstrict()
            @warn "Calculation of ỹ could be wrong due to insufficient tests."
        end
        z = y * inv(x)
        ỹ = tobits(z)[end]
        
        if mode == :compressed
            # [2/3; x]
            out[1] = ỹ ? 0x03 : 0x02
            copyto!(out, 2, _x, 1, nbytes_x)
            
        else # mode == :hybrid
            # [6/7; x; y]
            out[1] = ỹ ? 0x07 : 0x06
            copyto!(out, 2, _x, 1, nbytes_x)
            copyto!(out, nbytes_x + 2, _y, 1, nbytes_y)
        end
    else
        error("Unrecognized mode $mode")
    end
    
    return out
end


octet(p::AbstractPoint; mode::Symbol = :uncompressed) = _octet(gx(p), gy(p); mode)

iscompressable(::P) where P <: AbstractPoint = eq(P) <: Weierstrass 
