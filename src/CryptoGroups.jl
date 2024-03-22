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

# Some type piracy here
# Though this one can can be avoided if AbstractPoint would be defined in this module as well as BinaryField and primeField.

import .Specs: octet
using .Specs: point

octet(p::AbstractPoint; mode::Symbol = :uncompressed) = octet(value(gx(p)), value(gy(p)), spec(p); mode)
Base.convert(::Type{P}, po::Vector{UInt8}) where P <: AbstractPoint = P <| point(po, spec(P))
Base.rem(p::AbstractPoint, q::Integer) = rem(gx(p), q)


using .Specs: int2octet, bits2octet, BinaryBasis, octet2int
using .Fields: PrimeField, BinaryField

octet(x::BinaryField) = bits2octet(tobits(x))
Base.rem(x::BinaryField, q::Integer) = rem(octet(x) |> octet2int, q)

octet(x::PrimeField) = int2octet(value(x), bitlength(modulus(x)))
Base.rem(x::PrimeField, q::Integer) = rem(value(x), q)


Base.convert(group::Type{<:PGroup}, element::Vector{UInt8}) = group <| (element |> Specs.octet2int)
PGroup{S}(element::Vector{UInt8}) where S = convert(PGroup{S}, element)
octet(x::PGroup) = Specs.int2octet(value(x), bitlength(modulus(x)))
Base.rem(x::PGroup, q::Integer) = rem(value(x), q)

octet(g::ECGroup; mode = :uncompressed) = octet(g.x; mode)


export @bin_str

export spec, specialize, order, bitlength

import .Fields: Field, BinaryField, FP, F2GNB, F2PB, PrimeField
export Field, FP, F2GNB, F2PB 

import .Curves: AbstractPoint, ECPoint, AffinePoint, BinaryCurve, gx, gy, a, b, oncurve
export AbstractPoint, ECPoint, AffinePoint, BinaryCurve 

import .Specs: MODP, Koblitz, ECP, EC2N, PB, GNB, generator, Hash, PRG, RO, ROPRG, point, octet, curve
export MODP, Koblitz, ECP, EC2N, PB, GNB, generator, Hash, PRG, RO, ROPRG, Specs, curve

export ElGamal

end
