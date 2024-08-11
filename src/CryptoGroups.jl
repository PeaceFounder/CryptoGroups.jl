module CryptoGroups

global strict_mode::Bool = false

function set_strict_mode(mode::Bool)

    if strict_mode && !mode
        @warn "Strict mode for CryptoGroups is downgraded."
    end

    global strict_mode = mode

    return
end

isstrict() = strict_mode


include("utils.jl")

function spec end
function concretize_type end
function order end
function bitlength end
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

# A temporary declarations

import .Fields: octet, int2octet, octet2int, octet2bits, bits2octet

function point end
function hex2bits end

function concretize_type end

include("Specs/Specs.jl")
include("groups.jl")
include("spec.jl")

using .Fields: PrimeField, BinaryField

Base.rem(p::AbstractPoint, q::Integer) = rem(gx(p), q)
Base.rem(x::BinaryField, q::Integer) = rem(octet(x) |> octet2int, q)
Base.rem(x::PrimeField, q::Integer) = rem(value(x), q)
Base.rem(x::PGroup, q::Integer) = rem(value(x), q)
Base.rem(x::ECGroup, q::Integer) = rem(x.x, q)


point(po::String, spec::GroupSpec) = point(hex2bytes(po), spec)

function point(po::Vector{UInt8}, spec::GroupSpec)

    P = concretize_type(AffinePoint, spec)

    p = P(po)

    return point(p)
end

point(p::AffinePoint{<:Weierstrass}) = convert(Tuple{BigInt, BigInt}, p)
point(p::AffinePoint{<:BinaryCurve}) = convert(Tuple{BitVector, BitVector}, p)


export @bin_str

export spec, concretize_type #, order, bitlength

import .Fields: Field, BinaryField, FP, F2GNB, F2PB, PrimeField
#export Field, FP, F2GNB, F2PB 

import .Curves: AbstractPoint, ECPoint, AffinePoint, BinaryCurve, gx, gy, a, b, oncurve
export AbstractPoint, ECPoint, AffinePoint #, BinaryCurve 

import .Specs: generator
export generator, octet, order, modulus, value

end
