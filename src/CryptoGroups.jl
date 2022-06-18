module CryptoGroups


include("utils.jl")


function spec end
function specialize end
function order end
function bitlength end


<|(x::Type, y) = convert(x, y)


function tobits end

function value end
function modulus end
function cofactor end

function name end # Needed because of Curves.

function a end
function b end

# In most cases probably using sizeof directly is more appropriate
bitlength(::Type{T}) where T <: Integer = sizeof(T) * 8

function bitlength(p::Integer)

    bits = bitstring(p)
    start = findfirst(x -> x == '1', bits)
    N = length(bits) - start + 1

    return N
end


function bitlength(p::BigInt)

    bytes = int2bytes(p)
    bits = bitstring(bytes[end])
    start = findfirst(x -> x == '1', bits)
    N = length(bytes) * 8  - (start - 1)

    return N
end


include("Fields/Fields.jl")
include("Curves/Curves.jl")
include("Specs/Specs.jl")


include("groups.jl")
include("elgamal.jl")
include("spec.jl")


### Some type piracy here

import .Specs: octet
using .Specs: point

# Though this one can can be avoided if AbstractPoint would be defined in this module
octet(p::AbstractPoint; mode::Symbol = :uncompressed) = octet(value(gx(p)), value(gy(p)), spec(p); mode)
Base.convert(::Type{P}, po::Vector{UInt8}) where P <: AbstractPoint = P <| point(po, spec(P))

using .Specs: int2octet, bits2octet, BinaryBasis, octet2int, octet2bits
using .Fields: PrimeField, BinaryField

octet(x::BinaryField) = bits2octet(tobits(x))
octet(x::PrimeField) = int2octet(value(x), bitlength(modulus(x)))

# This for now is a temporary function

function ecp(; p::BigInt, n::BigInt, h::Int, a::Vector{UInt8}, b::Vector{UInt8}, G::Vector{UInt8})
    
    a_ = octet2int(a)
    b_ = octet2int(b)
    
    sp = ECP(p, n, a_, b_, nothing, nothing)

    #@infiltrate

    (Gx, Gy) = point(G, sp)

    return ECP(p, n, a_, b_, Gx, Gy)
end


function ec2n(basis::BinaryBasis; a::Vector{UInt8}, b::Vector{UInt8}, G::Vector{UInt8}, n::BigInt, h::Int)
    
    a_ = octet2bits(a, bitlength(basis))
    b_ = octet2bits(b, bitlength(basis))
    
    sp = EC2N(basis; a = a_, b = b_, n)

    (Gx, Gy) = point(G, sp)

    #(Gx, Gy) = reverse(Gx), reverse(Gy)

    return EC2N(basis, n, a_, b_, Gx, Gy)
end




export @bin_str

export spec, specialize, order, bitlength

import .Fields: Field, BinaryField, FP, F2GNB, F2PB, PrimeField
export Field, FP, F2GNB, F2PB 

import .Curves: AbstractPoint, ECPoint, AffinePoint, BinaryCurve, gx, gy, a, b, oncurve
export AbstractPoint, ECPoint, AffinePoint, BinaryCurve 

import .Specs: MODP, Koblitz, ECP, EC2N, PB, GNB, generator, Hash, PRG, RO, ROPRG, point, octet
export MODP, Koblitz, ECP, EC2N, PB, GNB, generator, Hash, PRG, RO, ROPRG

export ElGamal

end
