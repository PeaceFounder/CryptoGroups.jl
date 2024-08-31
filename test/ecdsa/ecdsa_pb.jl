using CryptoPRG: HashSpec
using CryptoGroups.Utils: @bin_str, @hex_str, octet2int
using CryptoGroups.Curves: gx, gy
using CryptoGroups.Specs: EC2N, PB
using CryptoGroups: concretize_type, ECPoint, generator, octet, order

# Elliptic Curve Domain Parameter Setup

basis = PB(hex"80000000 00000000 00000000 00000000 00000000 00000201", 191)

curve_spec = EC2N(basis;   ### Need to think about adding a proper methods
                  a = hex"2866537B 67675263 6A68F565 54E12640 276B649E F7526267",
                  b = hex"2E45EF57 1F00786F 67B0081B 9495A3D9 5462F5DE 0AA185EC",
                  cofactor = 2,
                  G = hex"04 36B3DAF8 A23206F9 C4F299D7 B21A9C36 9137F2C8 4AE1AA0D 765BE734 33B3F95E 332932E7 0EA245CA 2418EA0E F98018FB",
                  n = 1569275433846670190958947355803350458831205595451630533029,
)

P = concretize_type(ECPoint, curve_spec)
G = P(generator(curve_spec))

# Key Generation

d = 1275552191113212300012030439187146164646146646466749494799
Q = d*G


@test octet(Q) == hex"04 5DE37E75 6BD55D72 E3768CB3 96FFEB96 2614DEA4 CE28A2E7 55C0E0E0 2F5FB132 CAF416EF 85B229BB B8E13520 03125BA1"

# Signature Generation

M = "abc"
H = HashSpec("SHA1")
e = octet2int(H(M))

## Elliptic curve computaion

k = 1542725565216523985789236956265265265235675811949404040041
R = k*G

x = octet(gx(R)) # Instead I could have octet(gx(R))
y = octet(gy(R))

@test x == hex"438E5A11 FB55E4C6 5471DCD4 9E266142 A3BDF2BF 9D5772D5"
@test y == hex"2AD603A0 5BD1D177 649F9167 E6F475B7 E2FF590C 85AF15DA"

x̄ = octet2int(x)
ȳ = octet2int(y)

@test x̄ == 1656469817011541734314669640730254878828443186986697061077

n = order(P)
r = x̄ % n

@test r == 87194383164871543355722284926904419997237591535066528048

## Modular computation

s = invmod(k, n) * (e + d*r) % n

@test s == 308992691965804947361541664549085895292153777025772063598

# println("""The signature is:
#     r = $r
#     s = $s
# """)

# Signature verification

(r′, s′) = r, s

M′ = "abc"
H′ = HashSpec("sha256")
e′ = octet2int(H(M′))

@test 1 < r′ < n - 1
@test 1 < s′ < n - 1

c = invmod(s′, n)

@test c == 952933666850866331568782284754801289889992082635386177703


u₁ = e′*c % n
u₂ = r′*c % n

@test u₁ == 1248886407154707854022434516084062503301792374360994400066
@test u₂ == 527017380977534012168222466016199849611971141652753464154
    
W = u₁*G + u₂*Q 
x₁ = octet(gx(W))
y₁ = octet(gy(W))

@test x₁ == hex"438E5A11 FB55E4C6 5471DCD4 9E266142 A3BDF2BF 9D5772D5"
@test y₁ == hex"2AD603A0 5BD1D177 649F9167 E6F475B7 E2FF590C 85AF15DA"

x̄₁ = octet2int(x₁)
ν = x̄₁ % n

@test ν == r′
