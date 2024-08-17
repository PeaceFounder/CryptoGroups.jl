module Specs

# This module is concerned with loading exsiting or generating cryptiographic parameters

import ..CryptoGroups.Fields: modulus, bitlength, octet
import ..CryptoGroups.Curves: order, a, b, cofactor
import ..CryptoGroups: name #, point, octet
import ..CryptoGroups.Utils: octet2int, octet2bits, hex2bits, int2octet, bits2octet
import ..CryptoGroups: point

include("spec.jl")
include("field_specs.jl")

include("legacy.jl")

include("curve_constants.jl")
include("modp_constants.jl")

export ECP, EC2N, Koblitz, MODP

end
