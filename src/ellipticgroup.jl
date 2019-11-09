### Wrapping elliptic group from ECC
using ECC

struct EllipticGroup{T} <: CyclicGroup 
    G::T
    q
    t::Int
end

EllipticGroup{T}(sec,G::EllipticGroup{T}) where T <: S256Point = EllipticGroup(sec2point(sec),G.q,G.t)
EllipticGroup(sec,G::EllipticGroup) = typeof(G)(sec,G)

binary(G::EllipticGroup{T}) where T <: S256Point = point2sec(G.G)
value(G::EllipticGroup{T}) where T <: S256Point = G.G.𝑥.𝑛
security(G::EllipticGroup) = G.t
order(G::EllipticGroup) = G.q


function Scep256k1Group()
    G = ECC.G # G the scep256k1 generator point
    t = 256
    q = ECC.N
    EllipticGroup(G,q,t)
end

^(G::EllipticGroup,n) = EllipticGroup(n*G.G,G.q,G.t)

# When more elliptic groups appears then check for the elements should become important.
import Base.*
*(X::EllipticGroup,Y::EllipticGroup) = EllipticGroup(X.G + Y.G,X.q,X.t)

import Base.==
==(X::EllipticGroup,Y::EllipticGroup) = X.G==Y.G
