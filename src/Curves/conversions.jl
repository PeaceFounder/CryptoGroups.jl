# Implements point compression. Currently it is tied to AffinePoint but should be extended for generic case
using ..Fields: modulus, value, bitlength, int2octet, octet2int, octet2bits, bits2octet, tobits
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


function (::Type{P})(po::Vector{UInt8}) where P <: AbstractPoint

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


function _compressed_octet(x::Vector{UInt8}, y::Vector{UInt8}, ỹ::Bool)

    if ỹ == false
        return UInt8[2, x...]
    elseif ỹ == true
        return UInt8[3, x...]
    end

end


function _hybrid_octet(x::Vector{UInt8}, y::Vector{UInt8}, ỹ::Bool)

    if ỹ == false
        return UInt8[6, x..., y...]
    elseif ỹ == true
        return UInt8[7, x..., y...]
    end
    
end


function _uncompressed_octet(x::Vector{UInt8}, y::Vector{UInt8})
     return UInt8[4, x..., y...]
end


# This shall be considered internal as x, y can't be arbitrary!!!
function _octet(x::BigInt, y::BigInt, N::Int; mode::Symbol = :uncompressed) # N is bitlength(modulus(field()))

    _x = int2octet(x, N)
    _y = int2octet(y, N)

    if mode == :uncompressed

        return _uncompressed_octet(_x, _y)

    elseif mode in [:compressed, :hybrid]

        ỹ = mod(y, 2) % Bool
        
        if mode == :compressed
            
            return _compressed_octet(_x, _y, ỹ)
            
        elseif mode == :hybrid

            return _hybrid_octet(_x, _y, ỹ)

        end

    else
        error("Unrecognized mode $mode")
    end
    
end


_octet(x::F, y::F; mode::Symbol = :uncompressed) where F <: PrimeField = _octet(value(x), value(y), bitlength(modulus(F)); mode)


function _octet(x::F, y::F; mode::Symbol = :uncompressed) where F <: BinaryField

    _x = octet(x)
    _y = octet(y)

    if mode == :uncompressed

        return _uncompressed_octet(_x, _y)

    elseif mode in [:compressed, :hybrid]

        #x_ = F(reverse(x))
        #y_ = F(reverse(y))

        if isstrict()
            @warn "Calculation of ỹ could be wrong due to insufficient tests."
        end

        z = y * inv(x)

        ỹ = tobits(z)[end]

        if mode == :compressed
            
            return _compressed_octet(_x, _y, ỹ)
            
        elseif mode == :hybrid

            return _hybrid_octet(_x, _y, ỹ)

        end

    else
        error("Unrecognized mode $mode")
    end
end

octet(p::AbstractPoint; mode::Symbol = :uncompressed) = _octet(gx(p), gy(p); mode)

iscompressable(::P) where P <: AbstractPoint = eq(P) <: Weierstrass 
