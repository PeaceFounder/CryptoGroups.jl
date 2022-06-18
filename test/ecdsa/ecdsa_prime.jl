using Test
using CryptoGroups: @bin_str, @hex_str, ECP, specialize, ECPoint, <|, generator, octet, ecp, Hash, gx, gy, order, modinv
using CryptoGroups.Specs: octet2int


curve_spec = ecp(;

                 p = 6277101735386680763835789423207666416083908700390324961279,
                 n = 6277101735386680763835789423176059013767194773182842284081,
                 h = 1,
                 a = hex"FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE FFFFFFFF FFFFFFFC",
                 b = hex"64210519 E59C80E7 0FA7E9AB 72243049 FEB8DEEC C146B9B1",
                 G = hex"03 188DA80E B03090F6 7CBF20EB 43A18800 F4FF0AFD 82FF1012",
)


P = specialize(ECPoint, curve_spec)
G = P <| generator(curve_spec) 

# Key Generation

d = 651056770906015076056810763456358567190100156695615665659
Q = d*G


@test octet(Q; mode=:compressed) == hex"02 62B12D60 690CDCF3 30BABAB6 E69763B4 71F994DD 702D16A5"

# Signature Generation

M = "abc"
#H = Hash("sha256")
H = Hash("SHA1")
e = octet2int(H(M))

@test e == 968236873715988614170569073515315707566766479517

k = 6140507067065001063065065565667405560006161556565665656654

R = k*G

x₁ = octet(gx(R))
y₁ = octet(gy(R))

@test x₁ == hex"88505238 0FF147B7 34C330C4 3D39B2C4 A89F29B0 F749FEAD"
@test y₁ == hex"9CF9FA1C BEFEFB91 7747A3BB 29C072B9 289C2547 884FD835"


x̄₁ = octet2int(x₁)


@test x̄₁ == 3342403536405981729393488334694600415596881826869351677613

n = order(P)

r = x̄₁ % n

@test r == 3342403536405981729393488334694600415596881826869351677613

s = modinv(k, n) * (e + d*r) % n

@test s == 5735822328888155254683894997897571951568553642892029982342


println("""The signature is:
    r = $r
    s = $s
""")


# Signature verification

(r′, s′) = r, s

M′ = "abc"
H′ = Hash("sha256")
e′ = octet2int(H(M′))


@test 1 < r′ < n - 1
@test 1 < s′ < n - 1

c = modinv(s′, n)

@test c == 3250964404472526825130516490452346217749189704049629042861


u₁ = e′*c % n
u₂ = r′*c % n

@test u₁ == 2563697409189434185194736134579731015366492496392189760599
@test u₂ == 6266643813348617967186477710235785849136406323338782220568

W = u₁*G + u₂*Q 
x₁ = octet(gx(W))
y₁ = octet(gy(W))

@test x₁ == hex"88505238 0FF147B7 34C330C4 3D39B2C4 A89F29B0 F749FEAD"
@test y₁ == hex"9CF9FA1C BEFEFB91 7747A3BB 29C072B9 289C2547 884FD835"

x̄₁ = octet2int(x₁)

ν = x̄₁ % n

@test ν == r′

