#abstract type Spec end

module Specs


include("spec.jl")

include("field_specs.jl")

include("curve_constants.jl")
include("modp_constants.jl")


### For now keeping things simple
# Need to think about a more sustainable strategy
function spec(x::Symbol)
    if x == :P_192
        return Curve_P_192
    elseif x == :P_244
        return Curve_P_244
    elseif x == :P_256
        return Curve_P_256
    elseif x == :P_384
        return Curve_P_384
    elseif x == :P_521
        return Curve_P_521
    else
        error("$x not implemented")
    end
end


export ECP, EC2N, Koblitz, MODP, spec

end
