using Test
import CryptoGroups: PGroup, @PGroup, order, modulus, value, concretize_type, Specs, generator, octet

G = @PGroup{p = 23, q = 11}
q = order(G)
p = modulus(G)

@test isvalid(G(3)) == true
@test_throws ArgumentError G(11)

n = let 
    n = 0
    for i in 2:p-1
        isvalid(G(i; skip_validation=true)) && (n+=1)
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
import CryptoGroups: PGroup, concretize_type, value, order, ECGroup, Group, modulus, generator, Specs

using Test

### Testing prime groups
using Primes

function testgroup(g::G) where G <: Group
    @test g*g^2 == g^3
    @test (g^2)^2 == g^4
    # @test g^(order(G) + 1) == g # Done at the constructor

    @test one(g) * g == g
    @test g * one(g) == g

    @test convert(G, value(one(G)); allow_one=true) == one(G)
    @test convert(G, octet(one(G)); allow_one=true) == one(G)

end

G = @PGroup{p = 23, q = 11}
g = G(3)

@test value(g) == 3
@test order(G) == 11
@test modulus(G) == 23

testgroup(g)

### Testing ECGroup

let
    spec = Specs.Curve_P_256

    G = concretize_type(ECGroup, spec; name = :P_256)
    g = G(generator(spec))

    testgroup(g)
end

###  RFC standart group test with PGroup

let
    spec = Specs.MODP_1024
    G = concretize_type(PGroup, spec)
    g = G(generator(spec))

    testgroup(g)
end
### PGroup generation algorithms for 

using Random: MersenneTwister
rng = MersenneTwister(0)

# Broken
let 
modp_spec = Specs.sophie_germain_group(rng, 2, 100)
G = concretize_type(PGroup, modp_spec)
# g = G(generator(modp_spec))
# testgroup(g)
end

### Seems to have an issue with large numbers
let
modp_spec = Specs.dsa_standart_group(rng, 10, 10)
G = concretize_type(PGroup, modp_spec)
g = G(generator(modp_spec))

testgroup(g)
end

