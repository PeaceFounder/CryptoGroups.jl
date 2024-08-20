module CryptoGroups

import CryptoPRG: bitlength

global strict_mode::Bool = false

function set_strict_mode(mode::Bool)

    if strict_mode && !mode
        @warn "Strict mode for CryptoGroups is downgraded."
    end

    global strict_mode = mode

    return
end

isstrict() = strict_mode


function spec end
function concretize_type end
function order end
function name end # Needed because of Curves.
function value end

function iscompressable end

# Consider deprecating it in favour of value
# function point end

include("Utils.jl")
include("Fields/Fields.jl")
include("Curves/Curves.jl")
include("Specs/Specs.jl")
include("groups.jl")
include("spec.jl")
include("macros.jl")

using .Fields: PrimeField, BinaryField

Base.rem(p::AbstractPoint, q::Integer) = rem(gx(p), q)
Base.rem(x::BinaryField, q::Integer) = rem(octet(x) |> Utils.octet2int, q) # used in ec2n.jl test in CryptoSignatures
Base.rem(x::PrimeField, q::Integer) = rem(value(x), q)
Base.rem(x::PGroup, q::Integer) = rem(value(x), q)
Base.rem(x::ECGroup, q::Integer) = rem(x.x, q)

export spec, concretize_type #, order, bitlength

import .Specs: generator 
export generator, octet, order, value, @PGroup, @ECGroup, @ECPoint, Group

end
