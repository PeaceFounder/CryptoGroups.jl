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
@show value(G)
@show security(G)
@show order(G)
@show mod(G)
@show powermod(4,10,G)
@show inv(3,G)

testgroup(G)

### Testing EllitpicGroup

G = CryptoGroups.Scep256k1Group()
@show G^2
@show value(G)
@show security(G)
@show order(G)
@show mod(G)
@show powermod(4,10,G)
@show inv(3,G)

testgroup(G)

###  RFC standart

CryptoGroups.FirstOakleyGroup()
CryptoGroups.SecondOakleyGroup()
CryptoGroups.MODP160Group()
CryptoGroups.MODP224Group()
CryptoGroups.MODP256Group()

### Group generation algorithms

using Paillier
rng = Paillier.default_rng()

#G = CryptoGroups.SophieGermainGroup(rng,2,5)
#testgroup(G)

### Seems to be problem with large numbers
CryptoGroups.DSAStandartGroup(rng,10,10)
testgroup(G)
