### Wrapping elliptic group from ECC
using ECC
using BitConverter

struct EllipticGroup <: CyclicGroup # 
    G #::T
    q
    t::Int
end

#EllipticGroup{T}(sec::Vector{UInt8},G::EllipticGroup{T}) where T <: S256Point = EllipticGroup(sec2point(sec),G.q,G.t)
#EllipticGroup(sec::Vector{UInt8},G::EllipticGroup) = typeof(G)(sec,G)

EllipticGroup(n::BigInt, G::EllipticGroup) = EllipticGroup(sec2point(bytes(n)),G.q,G.t)

value(G::EllipticGroup) = to_big(point2sec(G.G))
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
