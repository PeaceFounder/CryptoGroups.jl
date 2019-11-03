### Wrapping elliptic group from ECC
using ECC

struct EllipticGroup <: AbstractGroup 
    G
    q
    t::Int
end

value(G::EllipticGroup) = G.G.𝑥.𝑛
security(G::EllipticGroup) = G.t
order(G::EllipticGroup) = G.q

function Scep256k1Group()
    G = ECC.G # G the scep256k1 generator point
    t = 256
    q = ECC.N
    EllipticGroup(G,q,t)
end

^(G::EllipticGroup,n) = EllipticGroup(n*G.G,G.q,G.t)
