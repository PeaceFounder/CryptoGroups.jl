abstract type Field end

(::Type{F})(x) where F <: Field = convert(F, x)

Base.:/(x::F, y::F) where F <: Field = x * inv(y)

Base.one(x::F) where F <: Field = one(F)
Base.zero(x::F) where F <: Field = zero(F)

square(x::Field) = x*x
Base.inv(x::Field) = x^(order(x) - 1)

Base.literal_pow(::typeof(^), x::F, ::Val{2}) where F <: Field = square(x)
Base.literal_pow(::typeof(^), x::F, ::Val{0}) where F <: Field = one(F)


abstract type BinaryField <: Field end

Base.:-(x::F, y::F) where F <: BinaryField = x + y

Base.convert(::Type{F}, x::BitVector) where F <: BinaryField = F(x)
#Base.convert(::Type{BitVector}, x::BinaryField) = 

function Base.convert(::Type{F}, x::Bool) where F <: BinaryField
    if x == false
        return zero(F)
    elseif x == true
        return one(F)
    end
end

Base.convert(::Type{F}, x::Integer) where F <: BinaryField = convert(F, Bool(x))

function Base.:^(g::F, a::Integer) where F <: BinaryField

    if a == 0
        return one(F)
    end

    e = tobits(a)

    
    r = length(e) - 1

    x = g

    for i in 2:r+1

        x = square(x)

        if e[i] == 1
            x = g*x
        end
    end

    return x
end


function (::Type{F})(x::Integer) where F <: BinaryField
    if x == 0
        return zero(F)
    elseif x==1
        return one(F)
    else
        error("Not supported")
    end
end

function Base.convert(::Type{F}, x) where F <: BinaryField 
    
    N = bitlength(F)

    bits = tobits(x)

    return convert(F, bits[end - N + 1:end])
end


Base.convert(::Type{F}, x::BinaryField) where F <: BinaryField = error("Can't be done blindly.")
Base.convert(::Type{F}, x::F) where F <: BinaryField = x

tobits(::Type{F}, x::BinaryField) where F <: BinaryField = convert(BitVector, x)

Base.show(io::IO, x::BinaryField) = print(io, join(i ? "1" : "0" for i in tobits(x)))

Base.isless(x::F, y::F) where F <: BinaryField = tobits(x) < tobits(y)

value(x::BinaryField) = tobits(x)

abstract type PrimeField <: Field end

Base.convert(::Type{F}, x::Integer) where F <: PrimeField = F(x)

function modulus end
function value end

Base.show(io::IO, x::PrimeField) = print(io, value(x))

modinv(s, q) = mod(gcdx(s, q)[2], q)

Base.inv(x::F) where F <: PrimeField = F(modinv(value(x), modulus(x)))

Base.zero(::Type{F}) where F <: PrimeField = F(0)
Base.one(::Type{F}) where F <: PrimeField = F(1)


Base.:+(x::F, y::F) where F <: PrimeField = F(mod(value(x) + value(y), modulus(F)))
Base.:-(x::F, y::F) where F <: PrimeField = x + F(modulus(F) - value(y))
Base.:-(x::F) where F <: PrimeField = F(modulus(F) - value(x))

Base.:*(x::F, y::F) where F <: PrimeField = F(mod(value(x) * value(y), modulus(F)))
Base.:*(x::Integer, y::F) where F <: PrimeField = F(x) * y
Base.:*(x::F, y::Integer) where F <: PrimeField = x * F(y)

Base.:^(x::F, n::Integer) where F <: PrimeField = F(powermod(value(x), n, modulus(F)))

Base.:(==)(x::F, y::F) where F <: PrimeField = value(x) == value(y)

Base.isless(x::F, y::F) where F <: PrimeField = value(x) < value(y)



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


octet(x::BinaryField) = bits2octet(tobits(x))
octet(x::PrimeField) = int2octet(value(x), bitlength(modulus(x)))

# Perhaps a convert method fits better here as the type is specific
(::Type{F})(x::Vector{UInt8}) where F <: PrimeField = F(octet2int(x))
(::Type{F})(x::Vector{UInt8}) where F <: BinaryField = F(octet2bits(x, bitlength(F)))
