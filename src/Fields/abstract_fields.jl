import ..CryptoGroups.Utils: int2octet!, octet2int, octet2bits, bits2octet

"""
    abstract type Field end

An abstract field interface. Requires `*`, `^`, `+`, `-`, `one`, `zero` and `order`, `octet`, `value`. Implements `inv`. The Field type
must supports a constructor for an element `e::Field` with properties:
```julia
e == typeof(e)(octet(e)) == typeof(e)(value(e))
```

"""
abstract type Field end

(::Type{F})(x) where F <: Field = convert(F, x)

Base.:/(x::F, y::F) where F <: Field = x * inv(y)

Base.one(x::F) where F <: Field = one(F)
Base.zero(x::F) where F <: Field = zero(F)

square(x::Field) = x*x
Base.inv(x::Field) = x^(order(x) - 1)

Base.literal_pow(::typeof(^), x::F, ::Val{2}) where F <: Field = square(x)
Base.literal_pow(::typeof(^), x::F, ::Val{0}) where F <: Field = one(F)

# Theese methods perhaps should be optionally enabled/disabled
Base.:*(k::Integer, x::F) where F <: Field = F(k) * x
Base.:*(x::Field, k::Integer) = k * x

Base.:/(x::F, k::Integer) where F <: Field = x / F(k)

Base.isnan(::Field) = false

"""
    abstract type BinaryField <: Field end

An interface for a a binary field. Requires constructor for `F(::BitVector)` and `convert(BitVector, ::BinaryField)`. Implements `octet` and coresponding constructor from octet. In addition implements a convinience constructor from integer by reinterpreting it first to bytes. 
"""
abstract type BinaryField <: Field end

Base.:-(x::F, y::F) where F <: BinaryField = x + y
Base.:-(x::BinaryField) = x 

Base.convert(::Type{F}, x::BitVector) where F <: BinaryField = F(x)

Base.eltype(::Type{F}) where F <: Field = F

function Base.convert(::Type{F}, x::Bool) where F <: BinaryField
    if x == false
        return zero(F)
    elseif x == true
        return one(F)
    end
end

function Base.convert(::Type{F}, x::Integer) where F <: BinaryField 

    if x == 0
        return zero(F)
    elseif x==1
        return one(F)
    else
        bytes = reinterpret(UInt8, [x])
        return F(reverse(bytes))
    end

end

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

(::Type{F})(x) where F <: BinaryField = convert(F, x)

function Base.convert(::Type{F}, x) where F <: BinaryField 
    
    N = bitlength(F)

    bits = tobits(x)

    return convert(F, bits[end - N + 1:end])
end


Base.convert(::Type{F}, x::BinaryField) where F <: BinaryField = error("Can't be done blindly.")
Base.convert(::Type{F}, x::F) where F <: BinaryField = x

tobits(::Type{F}, x::BinaryField) where F <: BinaryField = convert(BitVector, x)

function Base.show(io::IO, x::BinaryField)
    print(io, typeof(x))
    print(io, "(")
    print(io, "bin\"")
    print(io, join(i ? "1" : "0" for i in tobits(x)))
    print(io, "\")")
end

Base.isless(x::F, y::F) where F <: BinaryField = tobits(x) < tobits(y)

"""
    value(x::BinaryField)::BitVector

Converts a binary field elelement to a bitvector representation
"""
value(x::BinaryField) = tobits(x)

"""
    abstract type PrimeField <: Field end

An interface for a a modulus prime field. Requires constructor for `F(::BigInt)` and `convert(BigInt, ::PrimeField)` as well as `modulus`. It is recomended to subtype it for implementation of Mersenne primes.
"""
abstract type PrimeField <: Field end

Base.convert(::Type{F}, x::Integer) where F <: PrimeField = F(x)

"""
    modulus(::Union{F, Type{F}})::BigInt where F <: PrimeField

Gets a prime modulus of a 
"""
function modulus end

"""
    value(x::PrimeField)::BigInt

Converts a prime field element to an integer representation
"""
function value end

Base.show(io::IO, x::PrimeField) = print(io, value(x))

Base.inv(x::F) where F <: PrimeField = F(invmod(value(x), modulus(x)))

Base.zero(::Type{F}) where F <: PrimeField = F(0)
Base.one(::Type{F}) where F <: PrimeField = F(1)

"""
    +(x::F, y::F)::F where F <: Field

Adds two field elements.
"""
Base.:+(x::F, y::F) where F <: PrimeField = F(mod(value(x) + value(y), modulus(F)))

"""
    -(x::F)::F where F <: Field

Constructs a field element negative
"""
Base.:-(x::F, y::F) where F <: PrimeField = x + F(modulus(F) - value(y))
Base.:-(x::F) where F <: PrimeField = F(modulus(F) - value(x))

"""
    *(x::F, y::F)::F where F <: Field
    *(x::F, y::Integer)::F where F <: Field
    *(x::Integer, y::F)::F where F <: Field

Multiplies two field elements. If element is and integer it is first converted to the field element before used to multiply.
"""
Base.:*(x::F, y::F) where F <: PrimeField = F(mod(value(x) * value(y), modulus(F)))
Base.:*(x::Integer, y::F) where F <: PrimeField = F(x) * y
Base.:*(x::F, y::Integer) where F <: PrimeField = x * F(y)

"""
    ^(x::F, n::Integer) where F <: Field

Computes an exponential of a field element.
"""
Base.:^(x::F, n::Integer) where F <: PrimeField = F(powermod(value(x), n, modulus(F)))

Base.:(==)(x::F, y::F) where F <: PrimeField = value(x) == value(y)

Base.isless(x::F, y::F) where F <: PrimeField = value(x) < value(y)

"""
    octet(x::Field)

Returns a byte representation of a field element according to FIPS 186-4 standart.
"""
octet(x::BinaryField) = bits2octet(tobits(x))
#octet(x::PrimeField) = int2octet(value(x), bitlength(modulus(x)))
function octet(x::PrimeField)
    
    nbytes = cld(bitlength(x), 8)
    buffer = Vector{UInt8}(undef, nbytes)
    int2octet!(buffer, value(x), )

    return buffer
end


# Perhaps a convert method fits better here as the type is specific
(::Type{F})(x::Vector{UInt8}) where F <: PrimeField = F(octet2int(x))
(::Type{F})(x::Vector{UInt8}) where F <: BinaryField = F(octet2bits(x, bitlength(F)))


"""
    rem(x::BinaryField, q::T)::T where T <: Integer

Converts field value to an octet and then to integer from which remainder is computed. This is in accord with FIPS 186-4
standart for making ECDSA signatures over binary fields.
"""
Base.rem(x::BinaryField, q::Integer) = rem(octet(x) |> octet2int, q) # used in ec2n.jl test in CryptoSignatures

"""
    rem(x::PrimeField, q::T)::T where T <: Integer

Computes a remainder of the field integer value.
"""
Base.rem(x::PrimeField, q::Integer) = rem(value(x), q)
