module Curves

import ..CryptoGroups: <|, tobits, validate, a, b, cofactor, order, name

include("elliptic_curves.jl")
include("ecpoint.jl")

end
