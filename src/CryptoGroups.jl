module CryptoGroups

#using Infiltrator

include("utils.jl")

include("Specs/Specs.jl")

using .Specs: MODP, Koblitz, ECP, EC2N, Spec, PB, GNB
import .Specs: spec, order, a, b, bitlength, generator, modulus


include("elliptic_curves/abstract_fields.jl")
include("elliptic_curves/galois_fields.jl")
include("elliptic_curves/elliptic_curves.jl")
include("elliptic_curves/ecpoint.jl")
include("elliptic_curves/spec.jl")



include("groups/groups.jl")
include("groups/elgamal.jl")
include("groups/legacy.jl")
include("groups/primitives.jl") # May be better placed under `Specs.jl` to make a clear seperation from arithmetics and actual specification. 


export @bin_str

end
