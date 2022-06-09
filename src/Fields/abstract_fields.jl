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

Base.convert(::Type{F}, x::BitVector) where F <: BinaryField = F(x)

###
#function Base.convert(::Type{F}, x::Integer) where F <: BinaryField
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

