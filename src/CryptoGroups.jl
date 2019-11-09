module CryptoGroups

using Random: AbstractRNG

abstract type AbstractGroup end
abstract type CyclicGroup <: AbstractGroup end

import Base.^
^(G::AbstractGroup,n) = typeof(G)(G.G^n,G.q,G.t)

"""
The value of the group element.
"""
value(G::AbstractGroup) = error("Must be implemented by Group.")

"""
The number of bits one should use for the secret a and b in Diffie-Hellman key exchange.
"""
security(G::AbstractGroup) = error("Must be implemented by Group.")
order(G::AbstractGroup) = error("Must be implemented by Group.")
binary(G::AbstractGroup) = value(G) ### Default

include("utils.jl")
include("primegroup.jl")
include("rfc.jl")
include("ellipticgroup.jl")

export PrimeGroup, EllipticGroup, AbstractGroup, CyclicGroup, value, security, order, binary

end # module

