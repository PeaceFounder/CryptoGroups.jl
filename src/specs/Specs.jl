abstract type Spec end


specialize(::Type{ECGroup{P}}, spec::Spec; name = nothing) where P <: ECPoint = ECGroup{specialize(P, spec; name)}
specialize(::Type{ECGroup}, spec::Spec; name = nothing) = ECGroup{specialize(ECPoint, spec; name)}



function specialize(::Type{ECPoint{P}}, curve::Spec; name=nothing) where P <: AbstractPoint
    
    Q = specialize(P, curve)

    _order = order(curve)
    _cofactor = 1 # Need to update this one
    
    R = specialize(ECPoint{Q}, _order, _cofactor, name)

    return R
end

specialize(::Type{ECPoint}, spec::Spec; name = nothing) = specialize(ECPoint{AffinePoint}, spec; name)


##################### Macro for defining curve as also a group ###################

# macro def(constant_name, type, group_spec)

#     name_str = string(constant_name)

#     M = @__MODULE__

#     return esc(quote
#         const $constant_name = $M.specialize($type, $group_spec)

#         $M.name(::Type{$constant_name}) = $name_str

#         local ORDER = $group_spec.n
#         $M.order(::Type{$constant_name}) = ORDER

#         local GENERATOR = $constant_name($M.generator($group_spec)...)
#         $M.generator(::Type{$constant_name}) = GENERATOR
#     end)
# end

macro def(constant_name, type, group_spec)

    M = @__MODULE__

    return esc(quote
        const $constant_name = $M.specialize($type, $group_spec; name = $(QuoteNode(constant_name)))
    end)
end


include("curve_specs.jl")
include("modp_specs.jl")


