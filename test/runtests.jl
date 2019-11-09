using CryptoGroups
using Test

### Testing prime groups
using Primes

function testgroup(G)
    @test G*G^2 == G^3
    @test (G^2)^2 == G^4
    @test G^(order(G) + 1) == G
end

G = PrimeGroup(5,23,totient(23),4)
@show G^2
@show v = value(G)
@show security(G)
@show order(G)

@test PrimeGroup(binary(G^2),G)==G^2
@test typeof(G)(binary(G^2),G)==G^2

testgroup(G)

### Testing EllitpicGroup

G = CryptoGroups.Scep256k1Group()
@show G^2
@show value(G)
@show security(G)
@show order(G)

@test EllipticGroup(binary(G^2),G)==G^2
@test typeof(G)(binary(G^2),G)==G^2

testgroup(G)

###  RFC standart

CryptoGroups.FirstOakleyGroup()
CryptoGroups.SecondOakleyGroup()

G = CryptoGroups.MODP160Group()
testgroup(G)

G = CryptoGroups.MODP224Group()
testgroup(G)

G = CryptoGroups.MODP256Group()
#testgroup(G)

### Group generation algorithms

using Paillier
rng = Paillier.default_rng()

G = CryptoGroups.SophieGermainGroup(rng,2,100)
testgroup(G)

### Seems to be problem with large numbers
G = CryptoGroups.DSAStandartGroup(rng,10,10)
testgroup(G)
