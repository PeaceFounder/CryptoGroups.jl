# # Proof of Knowledge

# This Julia code example demonstrates the implementation of a non-interactive zero-knowledge proof of knowledge (NIZKPoK) for the discrete logarithm problem. The code showcases how to prove knowledge of an exponent $x$ for a given $y = g^x$ without revealing the secret $x$ itself, using group-agnostic methods.

# The implementation provides both simple proof by disclosure and a more secure non-interactive zero-knowledge proof, along with their corresponding verification functions. By leveraging Julia's multiple dispatch and parametric types, the code is designed to work with any abstract group `Group`, making it flexible and applicable to various cryptographic settings. For demonstration we use P-192 elliptic curve group, but the core functions are written to be compatible with any group provided by CryptoGroups.

using Test
using CryptoGroups
using CryptoGroups.Utils: int2octet
using CryptoPRG: HashSpec, Verificatum

verify(g::G, y::G, x::Integer) where G <: Group = g^x == y # proof by disclousure


function challenge(g::G, y::G, R::G) where G <: Group

    ## the encoding is deserializable as `octet` returns fixed length output that depends on unerlying group
    ## nevertheless it is recommended to use a proper canoncial encoding for this purpose which we shall skip
    prg = Verificatum.PRG(HashSpec("sha256"), [octet(g)..., octet(y)..., octet(R)...])

    return rand(prg, 2:order(G) - 1)
end

function prove(g::G, y::G, x::Integer) where G <: Group

    ## we can construct proof deterministically without relying on randomness source 
    prg = Verificatum.PRG(HashSpec("sha256"), [octet(y)..., int2octet(x)...]) 
    
    r = rand(prg, 2:order(G) - 1)

    R = g^r

    c = challenge(g, y, R)

    s = (r + c * x) % order(G) 

    return R, s
end

function verify(g::G, y::G, R::G, s::Integer) where G <: Group
    
    c = challenge(g, y, R)

    return g^s == R * y^c
end


g = @ECGroup{P_192}()

x = 21
y = g^x

@test verify(g, y, x)

proof = prove(g, y, x)
@test verify(g, y, proof...)
