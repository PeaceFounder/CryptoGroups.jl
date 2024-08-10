module Specs

# This module is concerned with loading exsiting or generating cryptiographic parameters

import ..CryptoGroups: order, a, b, cofactor, modulus, bitlength, name #, point, octet
import ..CryptoGroups: point, octet, octet2int, octet2bits, hex2bits, int2octet, bits2octet

include("spec.jl")
include("field_specs.jl")

include("legacy.jl")
#include("primitives.jl")

include("curve_constants.jl")
include("modp_constants.jl")

# Shall be refactored in Conversations instead
# include("conversions.jl")

export ECP, EC2N, Koblitz, MODP


end
