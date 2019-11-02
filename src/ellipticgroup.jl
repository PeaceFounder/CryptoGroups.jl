### Wrapping elliptic group from ECC
using ECC

struct EllipticGroup <: AbstractGroup 
    G
    t::Int
end

value(G::EllipticGroup) = G.G.𝑥.𝑛
security(G::EllipticGroup) = G.t

function Scep256k1Group()
    G = ECC.G # G the scep256k1 generator point
    t = 256
    EllipticGroup(G,t)
end

^(G::EllipticGroup,n) = EllipticGroup(n*G.G,G.t)
