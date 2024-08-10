using Test
import CryptoGroups: PGroup, order, modulus, value, specialize, Specs, generator

q = 11
p = 2*q + 1

G = PGroup(p, q; name=:G)

@test isvalid(G(3)) == true
@test isvalid(G(11)) == false

n = let 
    n = 0
    for i in 2:p-1
        isvalid(G(i)) && (n+=1)
    end
    n
end
@test n == q - 1

g = G(3)

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


@test one(g) * g == g
@test g * one(g) == g


################################ Legacy ##################################

import CryptoGroups
import CryptoGroups: PGroup, specialize, value, order, ECGroup, Group, modulus, generator, Specs

using Test

### Testing prime groups
using Primes

function testgroup(g::G) where G <: Group
    @test g*g^2 == g^3
    @test (g^2)^2 == g^4
    @test g^(order(G) + 1) == g

    @test one(g) * g == g
    @test g * one(g) == g
end

G = PGroup(23, 11)
g = G(3)

@test value(g) == 3
@test order(G) == 11
@test modulus(G) == 23

testgroup(g)

### Testing ECGroup

let
    spec = Specs.Curve_P_256

    G = specialize(ECGroup, spec; name = :P_256)
    g = G(generator(spec))

    testgroup(g)
end

###  RFC standart group test with PGroup

let
    spec = Specs.MODP_1024
    G = specialize(PGroup, spec)
    g = G(generator(spec))

    testgroup(g)
end
### PGroup generation algorithms for 

using Random: MersenneTwister
rng = MersenneTwister(0)

# Broken
let 
modp_spec = Specs.sophie_germain_group(rng, 2, 100)
G = specialize(PGroup, modp_spec)
g = G(generator(modp_spec))

#testgroup(g)
end

### Seems to have an issue with large numbers
let
modp_spec = Specs.dsa_standart_group(rng, 10, 10)
G = specialize(PGroup, modp_spec)
g = G(generator(modp_spec))

testgroup(g)
end

