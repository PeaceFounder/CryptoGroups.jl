module Fields

import ..CryptoGroups: bitlength, order, <|, tobits, value, modulus

include("abstract_fields.jl")
include("galois_fields.jl")

export Field, BinaryField, FP, F2PB, F2GNB

end
