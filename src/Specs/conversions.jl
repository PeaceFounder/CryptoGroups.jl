using ..CryptoGroups.Fields: F2PB, F2GNB, tobits
using CryptoUtils: sqrt_mod_prime


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

octet2int(x::String) =  octet2int(hex2bytes(x))

    
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


function decompress(x::BigInt, ỹ::Bool, spec::ECP)

    p = modulus(spec)

    α = mod(x^3 + a(spec) * x + b(spec), p)

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



function decompress(x::BitVector, ỹ::Bool, spec::EC2N)
    @warn "Not implemented, returning nothing..."
    return nothing
end


function point(x_octet::Vector{UInt8}, ỹ::Bool, spec::ECP)
    x = octet2int(x_octet)
    y = decompress(x, ỹ, spec)
    return (x, y)
end

function point(x::Vector{UInt8}, ỹ::Bool, spec::EC2N)
    error("Not implmented")
end


# I could also add validation that the point is on curve.
point(x::Vector{UInt8}, y::Vector{UInt8}, spec::ECP) = (octet2int(x), octet2int(y))

point(x::Vector{UInt8}, y::Vector{UInt8}, spec::EC2N) = (octet2bits(x, bitlength(spec)), octet2bits(y, bitlength(spec)))


function point(x_octet::Vector{UInt8}, y_octet::Vector{UInt8}, ỹ::Bool, spec::ECP)

    x, y = point(x_octet, y_octet, spec)

    y′ = decompress(x, ỹ, spec)
    
    if !(y == y′)
        @warn "Contradictory compression."
    end

    return (x, y)
end


function point(x_octet::Vector{UInt8}, y_octet::Vector{UInt8}, ỹ::Bool, spec::EC2N)

    x, y = point(x_octet, y_octet, spec)

    @warn "Sign bit is ignored."

    return (x, y)
end



function point(po::Vector{UInt8}, spec::Spec) 

    pc = po[1]

    if pc == 0
        return zero(spec)
    elseif pc == 2 
        ỹ = false
        x = po[2:end]
        return point(x, ỹ, spec)
    elseif pc == 3
        ỹ = true
        x = po[2:end]
        return point(x, ỹ, spec)
    elseif pc == 4 || pc == 6 || pc == 7
        
        l = div(length(po) - 1, 2)
        x = po[2:(l+1)]
        y = po[(l+2):(2l+1)]

        if pc == 4
            return point(x, y, spec)
        elseif pc == 6
            return point(x, y, false, spec)
        elseif pc == 7
            return point(x, y, true, spec)
        end
    else
        error("Wrong PC: $pc")
    end
end

point(po::String, spec::Spec) = point(hex2bytes(po), spec)


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




function octet(x::BigInt, y::BigInt, N::Int; mode::Symbol = :uncompressed)

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

octet(x::BigInt, y::BigInt, spec::ECP; mode::Symbol = :uncompressed) = octet(x, y, bitlength(modulus(spec)); mode)

_field(basis::PB) = F2PB(basis.f)
_field(basis::GNB) = F2GNB{basis.m, basis.T}


function octet(x::BitVector, y::BitVector, basis::BinaryBasis; mode::Symbol = :uncompressed)

    _x = bits2octet(x)
    _y = bits2octet(y)

    if mode == :uncompressed

        return _uncompressed_octet(_x, _y)

    elseif mode in [:compressed, :hybrid]

        F = _field(basis)

        # Reverse order is just a guess

        #x_ = F(x)
        #y_ = F(y)

        x_ = F(reverse(x))
        y_ = F(reverse(y))
        
        @warn "Calculation of ỹ could be wrong due to insufficient tests."

        z = y_ * inv(x_)

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

octet(x::BitVector, y::BitVector, spec::EC2N; mode::Symbol = :uncompressed) = octet(x, y, spec.basis; mode)
