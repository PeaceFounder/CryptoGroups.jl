module Curves

import ..CryptoGroups: order, name

include("elliptic_curves.jl")
include("ecpoint.jl")
include("conversions.jl")

export a, b, cofactor, gx, gy, order

end
