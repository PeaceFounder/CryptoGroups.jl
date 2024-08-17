using Test

import CryptoGroups
import CryptoGroups: generator, BinaryCurve, concretize_type, spec, Specs, @ECPoint
using CryptoGroups.Curves: order, oncurve, Weierstrass, ECPoint, AffinePoint, cofactor
using CryptoGroups.Fields: FP, F2PB, F2GNB



const FULL_TEST = false

for C in [:P_192, :P_224, :P_256, :P_384, :P_521]

    spec = getfield(Specs, Symbol("Curve_$C"))

    P = concretize_type(ECPoint{AffinePoint{Weierstrass, FP}}, spec)

    g = P(generator(spec))

    @test oncurve(g)
    @test oncurve(g*3)
    @test isvalid(g)

end


for C in [:B_163_PB, :B_233_PB, :B_283_PB, :B_409_PB, :B_571_PB, :K_163_PB, :K_233_PB, :K_283_PB, :K_409_PB, :K_571_PB]
    
    spec = getfield(Specs, Symbol("Curve_$C"))

    P = concretize_type(ECPoint{AffinePoint{BinaryCurve, F2PB}}, spec)

    g = P(generator(spec))

    @test oncurve(g)
    @test oncurve(g*3)
    
    if C == :K_163_PB # Too slow to test all 
        @test isvalid(g)
    end

end


for C in [:B_163_GNB, :B_233_GNB, :B_283_GNB, :B_409_GNB, :B_571_GNB, :K_163_GNB, :K_233_GNB, :K_283_GNB, :K_409_GNB, :K_571_GNB]

    spec = getfield(Specs, Symbol("Curve_$C"))

    P = concretize_type(ECPoint{AffinePoint{BinaryCurve, F2GNB}}, spec)

    g = P(generator(spec))

    @test oncurve(g)
    @test oncurve(g*3)

end


@test oncurve(@ECPoint{P_192}())
@test spec(@ECPoint{P_192}()) == Specs.Curve_P_192 # spec method constructs specification from a concrete type without lookup!
@test spec(@ECPoint{B_163}()) == Specs.Curve_B_163_GNB 
