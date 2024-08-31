module Fields

import ..CryptoGroups.Utils: tobits
import ..CryptoGroups: bitlength, order, value #, value, modulus

include("abstract_fields.jl")
include("galois_fields.jl")

export Field, FP, F2PB, @F2PB, F2GNB, value, modulus, octet

end
