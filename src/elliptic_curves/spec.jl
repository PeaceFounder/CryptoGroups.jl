# Need to introduce ECSpec to exclude MODP

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




function specialize(::Type{F2GNB}, N::Int)

    @assert div(N, 8) != 0 "Out of X9.62 spec"
    T = gn_basis_representation_rule(m)

    return F2GNB{N, T}
end

function F2GNB(x::BitVector)

    N = length(x)
    F = specialize(F2GNB, N)

    return F(x)
end


function specialize(::Type{AffinePoint{Weierstrass, F}}, curve::ECP) where F <: PrimeField

    (; p, a, b) = curve

    P = AffinePoint{specialize(Weierstrass, a, b), specialize(F, p)}
    
    return P
end

specialize(::Type{AffinePoint}, spec::ECP) = specialize(AffinePoint{Weierstrass, FP}, spec)

function specialize(::Type{F}, basis::PB) where F <: BinaryField
    (; f) = basis
    return specialize(F, f)
end


function specialize(::Type{F}, basis::GNB) where F <: BinaryField
    (; m, T) = basis
    return specialize(F, m, T)
end



specialize(::Type{BinaryCurve}, curve::EC2N) = specialize(BinaryCurve, a(curve), b(curve))

function specialize(::Type{AffinePoint{BinaryCurve, F}}, curve) where F <: BinaryField
    P = AffinePoint{specialize(BinaryCurve, curve), specialize(F, curve.basis)}
    return P
end

specialize(::Type{AffinePoint}, curve::EC2N{GNB}) = specialize(AffinePoint{BinaryCurve, F2GNB}, curve)
specialize(::Type{AffinePoint}, curve::EC2N{PB}) = specialize(AffinePoint{BinaryCurve, F2PB}, curve)



specialize(::Type{AffinePoint{BinaryCurve, F}}, curve::Koblitz) where F <: BinaryField = specialize(AffinePoint{BinaryCurve, F}, curve.bec)

function spec(::Type{P}; kwargs...) where P <: AbstractPoint

    n = hasmethod(order, Tuple{Type{P}}) ? order(P) : nothing
    h = hasmethod(cofactor, Tuple{Type{P}}) ? cofactor(P) : nothing

    EQ = eq(P)
    F = field(P)

    return spec(EQ, F; n, h, kwargs...)
end

spec(p::P) where P <: AbstractPoint = spec(P; Gx=value(gx(p)), Gy=value(gy(p)))


function spec(::Type{EQ}, ::Type{F}; n=nothing, h=nothing, Gx=nothing, Gy=nothing) where {EQ <: Weierstrass, F <: PrimeField}
    
    _a = a(EQ)
    _b = b(EQ)

    p = modulus(F)

    return ECP(p, n, _a, _b, Gx, Gy) # I will need to add H in the end here
end
