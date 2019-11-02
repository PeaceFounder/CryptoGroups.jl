module CryptoGroups

using Random: AbstractRNG

abstract type AbstractGroup end

import Base.^
^(G::AbstractGroup,n) = typeof(G)(G.G^n,G.t)

"""
The value of the group element.
"""
value(G::AbstractGroup) = error("Must be implemented by Group.")

"""
The number of bits one should use for the secret a and b in Diffie-Hellman key exchange.
"""
security(G::AbstractGroup) = error("Must be implemented by Group.")

include("utils.jl")
include("primegroup.jl")
include("rfc.jl")
include("ellipticgroup.jl")

export PrimeGroup, AbstractGroup, value, security

end # module

