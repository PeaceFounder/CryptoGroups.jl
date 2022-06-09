using ..CryptoGroups: static

struct ECPoint{P<:AbstractPoint, S} <: AbstractPoint # The same contract is satisfied thus a subtype
    p::P

    function ECPoint{P}(order::Integer, cofactor::Integer; name=nothing) where P <: AbstractPoint
        svars = static(; order, cofactor, name)
        return ECPoint{P, svars}
    end

    function ECPoint{P, S}(x::P) where {P <: AbstractPoint, S}
        
        @assert oncurve(x) "Point is not in curve"
        # A test with cofactor also here
        new{P, S}(x)
    end

    function ECPoint(p::P, order::Integer, cofactor::Integer; name=nothing) where P <: AbstractPoint
        EP = ECPoint{P}(order, cofactor; name)
        return EP(p)
    end
end

order(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = convert(Integer, S.order) 
cofactor(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = convert(Integer, S.cofactor) 

#modulus(::Type{ECPoint{P}}) where P <: AbstractPoint = modulus(P) # Only makes sense for point in prime fields!

name(::Type{ECPoint}) = nothing
name(::Type{ECPoint{P}}) where P <: AbstractPoint = nothing
name(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = isnothing(S.name) ? nothing : convert(Symbol, S.name)


eq(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = eq(P)
field(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = field(P)


### May need to do epoint seperatelly
function Base.show(io::IO, ::Type{P}) where P <: ECPoint
    if @isdefined P
        if name(P) == nothing
            Base._show_type(io, P)
        else
            print(io, name(P))
        end
    else
        print(io, "ECPoint")
    end
end


function Base.show(io::IO, p::P) where P <: ECPoint
    show(io, P)
    print(io, " <| (")
    show(io, gx(p))
    print(io, ", ")
    show(io, gy(p))
    print(io, ")")
end

function Base.display(::Type{P}) where P <: ECPoint
    show(P)
    ### I could be more precise on the Point
    ### Like ECPoint{AffinePoint{<:Weierstrass, <:F2GNB}, S}
    #if name(P) != nothing

    if !isnothing(name(P))

        print(" (alias for ECPoint{<:AffinePoint, ::static_ECPoint})") 
    end
end


function Base.:+(x::P, y::P) where P <: ECPoint 
    @assert x.p != y.p "AffinePoint's can't be equal. Multiply by 2 instead." 
    # Also assertion that x.p != -x.p
    return P(x.p + y.p)
end

Base.:*(x::P, n::Integer) where P <: ECPoint = P(x.p * n)
Base.:*(n::Integer, x::ECPoint) = x * n


Base.convert(::Type{ECPoint{P, S}}, x) where {P <: AbstractPoint, S} = ECPoint{P, S}(P <| x)
Base.convert(::Type{P}, x::P) where P <: ECPoint = x 


validate(p::P) where P <: ECPoint = order(P) * p.p == zero(p.p) 

oncurve(p::ECPoint) = oncurve(p.p)

Base.:(==)(x::ECPoint, y::ECPoint) = x.p == y.p


gx(p::ECPoint) = gx(p.p)
gy(p::ECPoint) = gy(p.p)


