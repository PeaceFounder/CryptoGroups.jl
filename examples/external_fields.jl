# # Field Subtyping

# This example demonstrates the flexibility and extensibility of the CryptoGroups package, showcasing how custom field implementations can be seamlessly integrated with existing elliptic curve operations. Specifically, it illustrates how to wrap and use `GaloisFields`'s binary extension field implementation within the `CryptoGroups` to do computations with elliptic curves.

# The modular approach of `CryptoGroups` allows users to integrate custom field implementations with elliptic curves easily. By defining a few key methods and ensuring proper subtyping, users can drop in their optimised field implementations - be it binary fields as shown here or Mersenne primes and test different optimisations with elliptic curves while allowing the use of existing implementations as reference. This design promotes the separation of concerns where one can focus on optimising fields, one new elliptic curve support or new arithmetic formulas like projective coordinates separately from each other. 

using Test
using CryptoGroups
import CryptoGroups: BinaryField, concretize_type, PB, bitlength
import GaloisFields: ExtensionField, GaloisField

struct GF₂{F <: ExtensionField} <: BinaryField
    x::F

    GF₂{F}(coeffs::NTuple{N, GaloisField(2)}) where {N, F <: ExtensionField} = new(F(coeffs))

    function GF₂{F}(xrev::BitVector) where F <: ExtensionField
        x = reverse(xrev)
        𝔽₂ = GaloisField(2)
        coeffs = Tuple(𝔽₂(i) for i in x) # May need a reverse
        return GF₂{F}(coeffs)
    end

    GF₂{F}(x::F) where F <: ExtensionField = new{F}(x)
    GF₂(x::F) where F <: ExtensionField = new{F}(x)
    
    function GF₂(f::BitVector)
        F, x = GaloisField(2, :x=>f)
        return GF₂{F}
    end

    GF₂(poly::Vector{Int}) = GF₂(BitVector(i in poly for i in 0:maximum(poly)))
end

bitlength(::Type{F}) where F <: ExtensionField = length(F.parameters[4]) - 1
bitlength(::Type{GF₂{F}}) where F <: ExtensionField = bitlength(F)

concretize_type(::Type{GF₂}, basis::PB) = GF₂(basis.f)

Base.convert(::Type{F}, a::BitVector) where F <: GF₂ = F(a)
Base.convert(::Type{BitVector}, a::GF₂) = reverse(BitVector(i.n for i in a.x.coeffs))

Base.:+(a::F, b::F) where F <: GF₂ = F(a.x + b.x)
Base.:*(a::F, b::F) where F <: GF₂ = F(a.x * b.x)

Base.inv(a::F) where F <: GF₂ = F(inv(a.x))

Base.zero(::Type{GF₂{F}}) where F <: ExtensionField = GF₂(zero(F))
Base.one(::Type{GF₂{F}}) where F <: ExtensionField = GF₂(one(F))


########################## This one we can test easally as follows ##################

import CryptoGroups.Curves: oncurve, order, AffinePoint, ECPoint, BinaryCurve
import CryptoGroups: generator
import CryptoGroups.Specs: Curve_B_163_PB, Curve_K_163_PB

let 
    B_163v3 = concretize_type(AffinePoint{BinaryCurve, GF₂}, Curve_B_163_PB)

    g = B_163v3(generator(Curve_B_163_PB))
    q = order(Curve_B_163_PB)

    @test oncurve(g)
    @test oncurve(g*3)
    @test g * (q + 1) == g
end

let
    K_163v3 = concretize_type(ECPoint{AffinePoint{BinaryCurve, GF₂}}, Curve_K_163_PB)

    g = K_163v3(generator(Curve_K_163_PB))
    q = order(Curve_K_163_PB)

    @test oncurve(g)
    @test oncurve(g*3)
    @test g * (q + 1) == g
end
