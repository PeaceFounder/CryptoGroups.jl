###################

using CryptoUtils: sqrt_mod_prime

# Field element to integer. Does not make sense for binary curves

<|(::Type{AffinePoint{EQ, F}}, args::Tuple{F, F}) where {EQ <: EllipticCurve, F <: Field} = AffinePoint{EQ, F}(args...)


function <|(::Type{P}, bytes::Tuple{Vector{UInt8}, Vector{UInt8}}) where P <: AbstractPoint
    xb, yb = bytes
    F = field(P)
    P <| (F <| xb, F <| yb)
end


function <|(::Type{P}, compressed_point::Tuple{Vector{UInt8}, UInt8}) where P <: AbstractPoint

    xb, ỹb = compressed_point
    
    F = field(P)
    
    x = F <| xb
    ỹ = Bool(ỹb)

    return P <| (x, ỹ)
end

#

function <|(::Type{AffinePoint{EQ, F}}, compressed_point::Tuple{BigInt, Bool}) where {EQ <: Weierstrass, F <: PrimeField}

    P = AffinePoint{EQ, F}

    x, ỹ = compressed_point

    p = modulus(F)

    α = mod(x^3 + a(EQ) * x + b(EQ), p)

    β = sqrt_mod_prime(α, p)

    # Rightmost bit of β. Perhaps this works:
    # Alternativelly I could do circshift
    if mod(β, 2) == ỹ
        y = β
    else
        y = p - β
    end

    return P <: (x, y)
end



function int2octet(x::Integer)

    hex = string(x, base=16)
    if mod(length(hex), 2) != 0
        hex = string("0", hex)
    end
    
    return hex2bytes(hex)
end


octet2int(x::Vector{UInt8}) = parse(BigInt, bytes2hex(x), base=16)



<|(::Type{Vector{UInt8}}, x::PrimeField) = int2octet(value(x))
<|(::Type{F}, x::Vector{UInt8}) where F <: PrimeField = F <| octet2int(x)
    

function octet2bits(x::Vector{UInt8})
    bv = BitVector(u << -i % Bool for u in x for i in 7:-1:0)
    return bv
end


octet2bits(x::Vector{UInt8}, N::Int) = octet2bits(x)[end - N + 1:end]


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


<|(::Type{Vector{UInt8}}, x::BinaryField) = bits2octet(tobits(x))
<|(::Type{F}, x::Vector{UInt8}) where F <: BinaryField = F <| octet2bits(x, bitlength(F))


function <|(::Type{P}, po::Vector{UInt8}) where P <: AbstractPoint

    pc = po[1]

    if pc == 0
        return zero(P)
    elseif pc == 2 
        ỹ = false
        x = po[2:end]
        return P <| (x, ỹ)
    elseif pc == 3
        ỹ = true
        x = po[2:end]
        return P <| (x, ỹ)
    elseif pc == 4 || pc == 6 || pc == 7
        l = div(length(po) - 1, 2)
        x = po[2:(l+1)]
        y = po[(l+2):(2l+1)]
        # Can write assertion for the rightmost bit of y here

        #@infiltrate
        return P <| (x, y)
    else
        error("Wrong PC: $pc")
    end

end


function <|(::Type{Vector{UInt8}}, p::AbstractPoint; option::Symbol=:uncompressed)

    if option == :uncompressed
        x = Vector{UInt8} <| gx(p)
        y = Vector{UInt8} <| gy(p)

        return UInt8[4, x..., y...]
    else
        error("To be implemented")
    end
end


function <|(x::Tuple{DataType, Symbol}, y) where T 
    return <|(x[1], y; option=x[2])
end

