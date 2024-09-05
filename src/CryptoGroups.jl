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

"""
    spec(::Union{G, Type{G}})::GroupSpec where G <: Group

Constructs a specification of a group type or instance. It's general intended use is for debugging purposes and serves as an inverse to `concretize_type` method. If called with a group instance the the value of it is used for a generator for the group specification otherwise it is left empty.

See also `concretize_type`
"""
function spec end


"""
    concretize_type(::Type{T}, args...) <: T where T

Constructs a concrete subtype of union all or abstract type `T`. The arguments `args` are used to set concrtete values for type parameters of `T`, effectivelly `concretizing` a generic abstract type. The resulting type then can be instanitated by a value. 

See also `spec`
"""
function concretize_type end

# HEllo world

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

import .Specs: generator
export spec, concretize_type #, order, bitlength
export octet, order, value, @PGroup, @ECGroup, @ECPoint, Group

end
