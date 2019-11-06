module CryptoGroups

using Random: AbstractRNG

abstract type AbstractGroup end
abstract type CyclicGroup <: AbstractGroup end

import Base.^
^(G::AbstractGroup,n) = typeof(G)(G.G^n,G.q,G.t)

import Base.mod
mod(G::CyclicGroup) = mod(value(G),order(G))
mod(n::Integer,G::CyclicGroup) = mod(n,order(G))

import Base.powermod
powermod(n::Integer,k::Integer,G::CyclicGroup) = powermod(n,k,order(G))

import Base.inv
inv(k::Integer,G::CyclicGroup) = powermod(k,order(G)-2,G)

"""
The value of the group element.
"""
value(G::AbstractGroup) = error("Must be implemented by Group.")

"""
The number of bits one should use for the secret a and b in Diffie-Hellman key exchange.
"""
security(G::AbstractGroup) = error("Must be implemented by Group.")
order(G::AbstractGroup) = error("Must be implemented by Group.")

include("utils.jl")
include("primegroup.jl")
include("rfc.jl")
include("ellipticgroup.jl")

export PrimeGroup, AbstractGroup, value, security, order

end # module

