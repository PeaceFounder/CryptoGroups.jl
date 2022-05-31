using Test
import CryptoGroups: PGroup, validate, order, Enc, Dec, modulus, value, specialize


q = 11
p = 2*q + 1


G = specialize(PGroup, p, q, :G)


@test validate(G(3)) == true
@test validate(G(11)) == false

n = let 
    n = 0
    for i in 2:p-1
        validate(G(i)) && (n+=1)
    end
    n
end
@test n == q - 1



g = G(3)

#q = 17
#g = PrimeGenerator(3, q)

@test g*g^2 == g^3
@test (g^2)^2 == g^4
@test g^(order(G) + 1) == g

h = g^7

@test h*h^2 == h^3
@test (h^2)^2 == h^4
@test h^(order(G) + 1) == h


@test inv(g)*g^2 == g
@test (g^7)^6 == g^(7*6) # This is only true for a cyclic group
@test g*g*g == g^3 # Checking multiplication
@test g^2/g == g



################################ Legacy ##################################

import CryptoGroups
import CryptoGroups: <|, PGroup, specialize, value, order, ECGroup, sophie_germain_group, dsa_standart_group, Group, modulus, generator, Specs

using Test

### Testing prime groups
using Primes

function testgroup(g::G) where G <: Group
    @test g*g^2 == g^3
    @test (g^2)^2 == g^4
    @test g^(order(G) + 1) == g
end

G = specialize(PGroup, 23, 11)
g = G(3)

@test value(g) == 3
@test order(G) == 11
@test modulus(G) == 23

testgroup(g)

### Testing ECGroup

let
    spec = Specs.Curve_P_256

    G = specialize(ECGroup, spec; name = :P_256)
    g = G <| generator(spec)

    testgroup(g)
end

###  RFC standart group test with PGroup

let
    spec = Specs.MODP_1024
    G = specialize(PGroup, spec)
    g = G <| generator(spec)

    testgroup(g)
end
### PGroup generation algorithms for 

using Random: MersenneTwister
rng = MersenneTwister(0)

# Broken
#g = sophie_germain_group(rng, 2, 100)
#testgroup(g)

### Seems to be problem with large numbers
g = dsa_standart_group(rng, 10, 10)
testgroup(g)


