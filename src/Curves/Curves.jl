module Curves

import ..CryptoGroups: order, name, value, concretize_type, iscompressable

include("ecpoint.jl")
include("elliptic_curves.jl")
include("conversions.jl")

export a, b, cofactor, gx, gy, order

end
