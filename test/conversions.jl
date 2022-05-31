# Point conversion routines
using Test
import CryptoGroups: bits2octet, octet2bits, spec, <|, AffinePoint, specialize, generator


let 
    x = UInt8[0, 129]
    @test bits2octet(octet2bits(x)) == x

    a = octet2bits(x)
    # First leftmostbits are padded to zero to make a full octet
    @test bits2octet(a[4:end]) == bits2octet(a)

    @test octet2bits(x, 13) == octet2bits(x)[4:end]
end

curve_spec = spec(:P_256)

P = specialize(AffinePoint, curve_spec)

p = P <| generator(curve_spec)

@test P <| Vector{UInt8} <| p == p

@test P <| (Vector{UInt8}, :uncompressed) <| p == p

