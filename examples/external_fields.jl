# This example shows how to use your own field implementation for using them to compute on 
# elliptic curves. For making things interesting we shall wrap GaloisFields implementeed binary extension field. 
using Test
using CryptoGroups
import CryptoGroups: BinaryField, tobits, frombits, specialize
import GaloisFields: ExtensionField, GaloisField

struct GF₂{F <: ExtensionField} <: BinaryField
    x::F
    GF₂{F}(coeffs::NTuple{N, GaloisField(2)}) where {N, F <: ExtensionField} = new(F(coeffs))

    function GF₂{F}(x::BitVector) where F <: ExtensionField
        𝔽₂ = GaloisField(2)
        coeffs = Tuple(𝔽₂(i) for i in x) # May need a reverse
        return GF₂{F}(coeffs)
    end

    GF₂{F}(x::F) where F <: ExtensionField = new{F}(x)
    GF₂(x::F) where F <: ExtensionField = new{F}(x)
end


function specialize(::Type{GF₂}, f::BitVector)
    F, x = GaloisField(2, :x=>f)
    return GF₂{F}
end

specialize(::Type{GF₂}, poly::Vector{Int}) = specialize(GF₂, BitVector(i in poly for i in 0:maximum(poly)))

tobits(a::GF₂) = reverse(BitVector(i.n for i in a.x.coeffs))
frombits(::Type{F}, a::BitVector) where F <: GF₂ = F(reverse(a))

Base.:+(a::F, b::F) where F <: GF₂ = F(a.x + b.x)
Base.:*(a::F, b::F) where F <: GF₂ = F(a.x * b.x)

Base.inv(a::F) where F <: GF₂ = F(inv(a.x))

Base.zero(::Type{GF₂{F}}) where F <: ExtensionField = GF₂(zero(F))
Base.one(::Type{GF₂{F}}) where F <: ExtensionField = GF₂(one(F))


########################## This one we can test easally as follows ##################

import CryptoGroups: Curve_B_163_PB, Curve_K_163_PB, BinaryCurve, generator, oncurve, order, @def, AffinePoint, <|, ECPoint, validate


let 
    B_163v3 = specialize(AffinePoint{BinaryCurve, GF₂}, Curve_B_163_PB)

    g = B_163v3 <| generator(Curve_B_163_PB)
    q = order(Curve_B_163_PB)

    @test oncurve(g)
    @test oncurve(g*3)
    @test g * q == zero(g)
end


@def K_163v3 ECPoint{AffinePoint{BinaryCurve, GF₂}} Curve_K_163_PB

let
    g = K_163v3 <| generator(Curve_K_163_PB)

    @test oncurve(g)
    @test oncurve(g*3)
    @test validate(g)
end
