using ..Fields: Field, PrimeField, BinaryField
import ..Fields

abstract type EllipticCurve end
abstract type AbstractPoint end

(::Type{P})(x) where P <: AbstractPoint = convert(P, x)


function (::Type{P})(x, y) where P <: AbstractPoint

    F = field(P)
    
    _x = F(x) # convert(F, x) could also be used
    _y = F(y)

    return P(_x, _y)
end


struct AffinePoint{E <: EllipticCurve, T <: Field} <: AbstractPoint
    x::T
    y::T

    AffinePoint{E, F}(x::F, y::F) where {E <: EllipticCurve, F <: Field} = new{E, F}(x, y)
    AffinePoint{E}(x::F, y::F) where {E <: EllipticCurve, F <: Field} = new{E, F}(x, y)
end

# AffinePoint{E, F}(x, y) where {E <: EllipticCurve, F <: Field} = AffinePoint{E, F}(convert(F, x), convert(F, y))

eq(::Type{AffinePoint{EQ, F}}) where {EQ <: EllipticCurve, F <: Field} = EQ
field(::Type{AffinePoint{EQ, F}}) where {EQ <: EllipticCurve, F <: Field} = F


### More through subtyping could be in here
Base.convert(::Type{P}, x::Tuple{BigInt, BigInt}) where P <: AffinePoint = P(x...)
Base.convert(::Type{P}, x::Tuple{BitVector, BitVector}) where P <: AffinePoint = P(x...)

Base.convert(::Type{Tuple{BigInt, BigInt}}, p::AbstractPoint) = (convert(BigInt, gx(p)), convert(BigInt, gy(p)))
Base.convert(::Type{Tuple{BitVector, BitVector}}, p::AbstractPoint) = (convert(BitVector, gx(p)), convert(BitVector, gy(p)))

gx(p::AffinePoint) = p.x
gy(p::AffinePoint) = p.y

Fields.modulus(::Type{AffinePoint{<:EllipticCurve, F}}) where F <: PrimeField = modulus(F)

Base.zero(::Type{AffinePoint{E, F}}) where {E <: EllipticCurve, F <: Field} = AffinePoint{E , F}(zero(F), zero(F))
Base.zero(x::P) where P <: AffinePoint = zero(P)

Base.:(==)(x::P, y::P) where P <: AffinePoint = x.x == y.x && x.y == y.y


Base.:-(u::P) where P <: AffinePoint{<:EllipticCurve, <: PrimeField} = P(u.x, -u.y)
Base.:-(u::P) where P <: AffinePoint{<:EllipticCurve, <: BinaryField} = P(u.x, u.x + u.y)

Base.:-(u::P, v::P) where P <: AffinePoint = u + (-v)

Base.isless(x::P, y::P) where P <: AbstractPoint = gx(x) == gx(y) ? gx(x) < gx(y) : gy(x) < gy(y)


function Base.:*(P::AffinePoint, k::Integer)
    
    h = tobits(3*k) 
    r = length(h) - 1
    e = tobits(k, length(h))

    R = P
    
    for i in 2:r #+1
        
        R = double(R) #R + R
        
        if h[i] == 1 && e[i] == 0
            R = R + P
        elseif h[i] == 0 && e[i] == 1

            if R == -P
                R = double(R)
            else
                R = R - P
            end

        end
    end

    return R
end

Base.:*(k::Integer, P::AffinePoint) = P * k


function Base.isvalid(::Type{P}) where P <: AffinePoint
    
    g = generator(x)

    if oncurve(g) == false
        return false
    end

    if !(g * order(P) == zero(P))
        return false
    end
    
    return true
end

Base.isvalid(x::AffinePoint) = oncurve(x)

### Definition of some elliptic curves and coresponding operations

struct Weierstrass{a, b} <: EllipticCurve end # May assume that a, b are 


a(::Type{Weierstrass{A, B}}) where {A, B} = A
b(::Type{Weierstrass{A, B}}) where {A, B} = B

a(::Type{AffinePoint{W, F}}) where {W <: Weierstrass, F <: Field} = convert(F, a(W))
b(::Type{AffinePoint{W, F}}) where {W <: Weierstrass, F <: Field} = convert(F, b(W))


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

    λ = (3 * x^2 + a(P))/ (2 * y) 
    
    x₃ = λ^2 - 2x
    y₃ = λ*(x - x₃) - y

    return P(x₃, y₃)
end

function oncurve(u::P) where P <: AffinePoint{<:Weierstrass, <:Field}
    (; x, y) = u
    return y^2 == x^3 + a(P)*x + b(P)
end
 

struct BinaryCurve{a, b} <: EllipticCurve end


a(::Type{BinaryCurve{A, B}}) where {A, B} = A
b(::Type{BinaryCurve{A, B}}) where {A, B} = B

a(::Type{AffinePoint{W, F}}) where {W <: BinaryCurve, F <: Field} = convert(F, a(W))
b(::Type{AffinePoint{W, F}}) where {W <: BinaryCurve, F <: Field} = convert(F, b(W))


function Base.:+(u::AffinePoint{E, F}, v::AffinePoint{E, F}) where {E<:BinaryCurve, F<:BinaryField}
        
    if u.x == v.x && u.y + v.y + u.x == zero(F)
        return zero(AffinePoint{E, F})
    end

    @assert u.x != v.x 

    _a = a(AffinePoint{E, F})

    λ = (v.y + u.y)/(v.x + u.x)

    x₃ = λ^2 + λ + u.x + v.x + _a
    y₃ = λ * (u.x + x₃) + x₃ + u.y

    return AffinePoint{E, F}(x₃, y₃)
end


function double(u::P) where P <: AffinePoint{<:BinaryCurve, <:BinaryField}
    
    (; x, y) = u
    
    λ = x + y/x
    
    x₃ = λ^2 + λ + a(P)
    y₃ = x^2 + (λ + one(λ))*x₃

    return P(x₃, y₃)
end


function oncurve(u::P) where P <: AffinePoint{<:BinaryCurve, <:Field}
    (; x, y) = u

    return y^2 + x*y == x^3 + a(P)*x^2 + b(P)
end
