module Specs

import ..CryptoGroups: order, a, b, cofactor, modulus, bitlength

include("spec.jl")
include("field_specs.jl")

include("legacy.jl")
include("primitives.jl")

include("curve_constants.jl")
include("modp_constants.jl")

export ECP, EC2N, Koblitz, MODP

end
