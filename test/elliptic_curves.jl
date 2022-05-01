using Test
using CryptoGroups
import CryptoGroups: F2PB, F2GNB, FP, Reducer, Weierstrass, AffinePoint, double, Reducer, BinaryCurve, oncurve, specialize

eq = Weierstrass{1, 1}

x1 = FP{23}(3)
y1 = FP{23}(10)

p1 = AffinePoint{eq}(x1, y1)

x2 = FP{23}(9)
y2 = FP{23}(7)

p2 = AffinePoint{eq}(x2, y2)

@test p1 + p2 == AffinePoint{eq, FP{23}}(17, 20)
@test double(p1) == AffinePoint{eq, FP{23}}(7, 12)
@test oncurve(p1) == true
@test oncurve(p2) == true

f = reverse(bin"10011")
R = Reducer(f)

α = F2PB{R}(bin"0100")

a = α^4
b = one(α)

#eq = BEQ(a, b)
eq = specialize(BinaryCurve, a, b)

p1 = AffinePoint{eq}(α^6, α^8)
p2 = AffinePoint{eq}(α^3, α^13)

@test p1 + p2 == AffinePoint{eq}(one(α), α^13)
@test double(p1) == AffinePoint{eq}(α^10, α^8)

@test oncurve(p1) == true
@test oncurve(p2) == true

α = F2GNB{4, 1}(bin"1100")

#eq = BEQ(zero(α), α^3)

eq = specialize(BinaryCurve, zero(α), α^3)

#error("here")

G = AffinePoint{eq}(α^3, α^5)


@test double(G) == AffinePoint{eq}(α^4, α^3)
@test double(G) + G == AffinePoint{eq}(α^13, α^2)

@test oncurve(G) == true
@test oncurve(G*3) == true



### Experimenting with 𝔽₂ field in Weierstrass equation. The doubling operation though does not work and thus is type constrained.

α = F2GNB{4, 1}(bin"1100")

#eq = WeierstrassEQ(one(α), one(α))
eq = specialize(Weierstrass, one(α), one(α))


p1 = AffinePoint{eq}(α^5, α^10)
p2 = AffinePoint{eq}(α^14, α^5)

@test oncurve(p1) == true
@test oncurve(p2) == true

@test oncurve(p1 + p2) == true

# Testing point multiplication
eq = Weierstrass{1, 1}

p = AffinePoint{eq, FP{23}}(3, 10)

@test oncurve(3*p) == true
@test double(p) == 2p
@test double(p) + p == 3p

@test double(p) + p + double(p) == 5p
@test double(double(p)) + p == 5p



############################### Real curve examples #####################

##################### P_192 ####
let

    p = 6277101735386680763835789423207666416083908700390324961279
    n = 6277101735386680763835789423176059013767194773182842284081

    a = p - 3 
    b = parse(BigInt, "64210519 e59c80e7 0fa7e9ab 72243049 feb8deec c146b9b1", base=16)

    Gx = parse(BigInt, "188da80e b03090f6 7cbf20eb 43a18800 f4ff0afd 82ff1012", base=16)
    Gy = parse(BigInt, "07192b95 ffc8da78 631011ed 6b24cdd5 73f977a1 1e794811", base=16)

    #point = AffinePoint{WeierstrassEQ(a, b), (FP/p)}(Gx, Gy)
    point = AffinePoint{specialize(Weierstrass, a, b), specialize(FP, p)}(Gx, Gy)


    @test oncurve(point) == true
    @test point * n == zero(point)

end

######################  K_163 #################

import CryptoGroups: hex2bits, <|

let

    a = 1
    b = 1

    ########### Polynomial basis

    #R = Reducer([163, 7, 6, 3, 0])
    #F = F2PB{R}
    
    F = specialize(F2PB, [163, 7, 6, 3, 0])

    EQ = specialize(BinaryCurve, F <| a, F <| b)

    Gx = "2 fe13c053 7bbc11ac aa07d793 de4e6d5e 5c94eee8" 
    Gy = "2 89070fb0 5d38ff58 321f2e80 0536d538 ccdaa3d9"

    point = AffinePoint{EQ, F}(Gx, Gy)
    @test oncurve(point)


    ######### Normal Basis

    T = 4 

    Gx = "0 5679b353 caa46825 fea2d371 3ba450da 0c2a4541" 
    Gy = "2 35b7c671 00506899 06bac3d9 dec76a83 5591edb2"

    F = F2GNB{163, 4}
    eq = specialize(BinaryCurve, F <| a, F <| b)

    point = AffinePoint{eq, F}(Gx, Gy)

    @test oncurve(point)
    @test oncurve(point*3) 

end
