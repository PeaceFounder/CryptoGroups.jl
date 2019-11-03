using CryptoGroups

### Testing prime groups

G = PrimeGroup(5,23,0,4)
@show G^2
@show value(G)
@show security(G)

### Testing EllitpicGroup

G = CryptoGroups.Scep256k1Group()
@show G^2
@show value(G)
@show security(G)

###  RFC standart

CryptoGroups.FirstOakleyGroup()
CryptoGroups.SecondOakleyGroup()
CryptoGroups.MODP160Group()
CryptoGroups.MODP224Group()
CryptoGroups.MODP256Group()

### Group generation algorithms

using Paillier
rng = Paillier.default_rng()

CryptoGroups.SophieGermainGroup(rng,2,256)

CryptoGroups.DSAStandartGroup(rng,100,100)
