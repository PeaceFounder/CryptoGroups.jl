module CryptoGroups

#using Infiltrator

include("utils.jl")

include("elliptic_curves/abstract_fields.jl")
include("elliptic_curves/galois_fields.jl")
include("elliptic_curves/elliptic_curves.jl")
include("elliptic_curves/basis_specs.jl")
include("elliptic_curves/ecpoint.jl")


include("groups/groups.jl")
include("groups/elgamal.jl")
include("groups/primitives.jl")
include("groups/legacy.jl")


include("specs/Specs.jl")


export @bin_str

end
