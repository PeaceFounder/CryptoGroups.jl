abstract type Field end

Base.:/(x::F, y::F) where F <: Field = x * inv(y)

Base.one(x::F) where F <: Field = one(F)
Base.zero(x::F) where F <: Field = zero(F)

square(x::Field) = x*x
Base.inv(x::Field) = x^(order(x) - 1)

Base.literal_pow(::typeof(^), x::F, ::Val{2}) where F <: Field = square(x)
Base.literal_pow(::typeof(^), x::F, ::Val{0}) where F <: Field = one(F)

abstract type BinaryField <: Field end

Base.:-(x::F, y::F) where F <: BinaryField = x + y


<|(::Type{F}, x) where F <: BinaryField = frombits(F, x)

function <|(::Type{F}, x::Integer) where F <: BinaryField
    if x == 0
        return zero(F)
    elseif x == 1
        return one(F)
    else
        error("Conversion of $x not possible")
    end
end


function Base.:^(g::F, a::Integer) where F <: BinaryField

    if a == 0
        return one(F)
    end

    e = tobin(a)
    
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


frombits(::Type{F}, gx::String) where F <: BinaryField = frombits(F, _hex2bytes(gx))
frombits(::Type{F}, gx::Vector{UInt8}) where F <: BinaryField = frombits(F, bytes2bits(gx))
frombits(::Type{F}, gx) where F <: BinaryField = frombits(F, convert(BitVector, gx))

frombits(::Type{F}, gx::BitVector) where F <: BinaryField = error("Must be implemented by field")


Base.show(io::IO, x::BinaryField) = print(io, join(i ? "1" : "0" for i in tobits(x)))

# Note if it becomes a bottleneck specialized methods can be added
Base.isless(x::F, y::F) where F <: BinaryField = tobits(x) < tobits(y)


abstract type PrimeField <: Field end


<|(::Type{F}, x) where F <: PrimeField = F(x)
<|(::Type{F}, x::StaticBigInt) where F <: PrimeField = F(BigInt(x))


function modulus end
function value end

Base.show(io::IO, x::PrimeField) = print(io, value(x))

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

