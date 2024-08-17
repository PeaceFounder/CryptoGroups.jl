module Curves

import ..CryptoGroups: order, name, value, concretize_type

include("elliptic_curves.jl")
include("ecpoint.jl")
include("conversions.jl")

export a, b, cofactor, gx, gy, order

end
