abstract type EllipticCurve end
abstract type AbstractPoint end



struct AffinePoint{E <: EllipticCurve, T <: Field} <: AbstractPoint
    x::T
    y::T

    AffinePoint{E, F}(x::F, y::F) where {E <: EllipticCurve, F <: Field} = new{E, F}(x, y)
    AffinePoint{E}(x::F, y::F) where {E <: EllipticCurve, F <: Field} = new{E, F}(x, y)
end


AffinePoint{E, F}(x, y) where {E <: EllipticCurve, F <: Field} = AffinePoint{E, F}(F <| x, F <| y)


eq(::Type{AffinePoint{EQ, F}}) where {EQ <: EllipticCurve, F <: Field} = EQ
field(::Type{AffinePoint{EQ, F}}) where {EQ <: EllipticCurve, F <: Field} = F


<|(::Type{P}, x::Tuple{BigInt, BigInt}) where P <: AffinePoint = P(x...)
<|(::Type{P}, x::Tuple{BitVector, BitVector}) where P <: AffinePoint = P(x...)

gx(p::AffinePoint) = p.x
gy(p::AffinePoint) = p.y

modulus(::Type{AffinePoint{<:EllipticCurve, F}}) where F <: PrimeField = modulus(F)

Base.zero(::Type{AffinePoint{E, F}}) where {E <: EllipticCurve, F <: Field} = AffinePoint{E , F}(zero(F), zero(F))
Base.zero(x::P) where P <: AffinePoint = zero(P)

Base.:(==)(x::P, y::P) where P <: AffinePoint = x.x == y.x && x.y == y.y


Base.:-(u::P) where P <: AffinePoint{<:EllipticCurve, <: PrimeField} = P(u.x, -u.y)
Base.:-(u::P) where P <: AffinePoint{<:EllipticCurve, <: BinaryField} = P(u.x, u.x + u.y)

Base.:-(u::P, v::P) where P <: AffinePoint = u + (-v)


function Base.:*(P::AffinePoint, k::Integer)
    
    h = tobin(3*k) 
    r = length(h) - 1
    e = tobin(k, length(h))

    R = P
    
    for i in 2:r #+1
        
        R = double(R) #R + R
        
        if h[i] == 1 && e[i] == 0
            R = R + P
        elseif h[i] == 0 && e[i] == 1
            R = R - P
        end
    end

    return R
end

Base.:*(k::Integer, P::AffinePoint) = P * k


function validate(::Type{P}) where P <: AffinePoint
    
    g = generator(x)

    if oncurve(g) == false
        return false
    end

    if !(g * order(P) == zero(P))
        return false
    end
    
    return true
end

validate(x::AffinePoint) = oncurve(x)


### Definition of some elliptic curves and coresponding operations

struct Weierstrass{a, b} <: EllipticCurve end # May assume that a, b are 

static_int(x::Integer) = x
static_int(x::BigInt) = StaticBigInt(x)

specialize(::Type{Weierstrass}, a::Integer, b::Integer) = Weierstrass{static_int(a), static_int(b)}

specialize(::Type{Weierstrass}, a::BitVector, b::BitVector) = Weierstrass{StaticBitVector(a), StaticBitVector(b)}
specialize(::Type{Weierstrass}, a::F, b::F) where F <: BinaryField = specialize(Weierstrass, tobits(a), tobits(b))


a(::Type{Weierstrass{A, B}}) where {A, B} = A
b(::Type{Weierstrass{A, B}}) where {A, B} = B

_a(::Type{AffinePoint{W, F}}) where {W <: Weierstrass, F <: Field} = F <| a(W)
_b(::Type{AffinePoint{W, F}}) where {W <: Weierstrass, F <: Field} = F <| b(W)


function Base.:+(u::AffinePoint{E, F}, v::AffinePoint{E, F}) where {E <: Weierstrass, F <: Field}
        
    if u.x == v.x && u.y + v.y == zero(F) # I could put this case into ECPoint as assertion
        return zero(AffinePoint{E, F})
    end

    @assert u.x != v.x 
    
    λ = (v.y - u.y)/(v.x - u.x)

    x₃ = λ^2 - u.x - v.x
    y₃ = λ * (u.x - x₃) - u.y

    return AffinePoint{E, F}(x₃, y₃)
end


function double(u::P) where P <: AffinePoint{<:Weierstrass, <:PrimeField}
    
    (; x, y) = u

    a = _a(P)

    λ = (3 * x^2 + a)/ (2 * y) 
    
    x₃ = λ^2 - 2x
    y₃ = λ*(x - x₃) - y

    return P(x₃, y₃)
end

function oncurve(u::P) where P <: AffinePoint{<:Weierstrass, <:Field}
    (; x, y) = u

    a = _a(P)
    b = _b(P)

    return y^2 == x^3 + a*x + b
end
 

struct BinaryCurve{a, b} <: EllipticCurve end

specialize(::Type{BinaryCurve}, a::BitVector, b::BitVector) = BinaryCurve{StaticBitVector(a), StaticBitVector(b)}
specialize(::Type{BinaryCurve}, a::F, b::F) where F <: BinaryField = specialize(BinaryCurve, tobits(a), tobits(b)) 


a(::Type{BinaryCurve{A, B}}) where {A, B} = A
b(::Type{BinaryCurve{A, B}}) where {A, B} = B

_a(::Type{AffinePoint{W, F}}) where {W <: BinaryCurve, F <: Field} = F <| a(W)
_b(::Type{AffinePoint{W, F}}) where {W <: BinaryCurve, F <: Field} = F <| b(W)


function Base.:+(u::AffinePoint{E, F}, v::AffinePoint{E, F}) where {E<:BinaryCurve, F<:BinaryField}
        
    if u.x == v.x && u.y + v.y + u.x == zero(F)
        return zero(AffinePoint{E, F})
    end

    @assert u.x != v.x 

    a = _a(AffinePoint{E, F})

    λ = (v.y + u.y)/(v.x + u.x)

    x₃ = λ^2 + λ + u.x + v.x + a
    y₃ = λ * (u.x + x₃) + x₃ + u.y

    return AffinePoint{E, F}(x₃, y₃)
end


function double(u::P) where P <: AffinePoint{<:BinaryCurve, <:BinaryField}
    
    (; x, y) = u

    a = _a(P)
    
    λ = u.x + u.y/u.x
    
    x₃ = λ^2 + λ + a
    y₃ = u.x^2 + (λ + one(λ))*x₃

    return P(x₃, y₃)
end


function oncurve(u::P) where P <: AffinePoint{<:BinaryCurve, <:Field}
    (; x, y) = u

    a = _a(P)
    b = _b(P)

    return y^2 + x*y == x^3 + a*x^2 + b
end



###################



using CryptoUtils: sqrt_mod_prime

### Here may be a good place for setting up conversion

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

# binary_order is abad name. Better would be to use already existing name nbits; Also applicable to field element
<|(::Type{F}, x::Vector{UInt8}) where F <: BinaryField = F <| octet2bits(x, binary_order(F))


#<|(::Type{F}, x::String) where F <: Field = F <| _hex2bytes(x) # 



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


### Should pass :hybrid as optional argument
### I could add similar arguments for integers
# F <|(:big) x or F <|(:little) x 
# 
#Vector{UInt8} <|(:hybrid) po
#Vector{UInt8} <|(:compressed) po
#Vector{UInt8} <|(:uncompressed) po

#Vector{UInt8} <| po % :uncompressed

###

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


# function f(x::Tuple{Type{T}, Symbol}, y) where T
#     return T(y)
# end

# f((BigInt, :hello), 4)


# function g(x::Type{T}, y) where T
#     return T(y)
# end

# g(BigInt, 4)


# function g(x::DataType, y) where T
#     return x[1](y)
# end
