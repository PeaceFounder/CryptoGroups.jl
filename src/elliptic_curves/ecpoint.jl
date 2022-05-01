# Static fields for ECPoint
struct static_ECPoint{N}
    order::StaticBigInt{N}
    cofactor::Int
    name::UInt128

    static_ECPoint(order::StaticBigInt{N}, cofactor::Int, name::Symbol) where N = new{N}(order, cofactor, string2uint(string(name)))
    static_ECPoint(order::Integer, cofactor::Integer, name::Symbol) = static_ECPoint(StaticBigInt(order), Int(cofactor), name)
    static_ECPoint(order::Integer, cofactor::Integer, ::Nothing) = static_ECPoint(StaticBigInt(order), Int(cofactor), Symbol(""))
end

struct ECPoint{P<:AbstractPoint, S} <: AbstractPoint # The same contract is satisfied thus a subtype
    p::P

    function ECPoint{P, S}(x::P) where {P <: AbstractPoint, S}
        
        @assert oncurve(x) "Point is not in curve"
        # A test with cofactor also here

        new{P, S}(x)
    end
end

specialize(::Type{ECPoint{P}}, order::Integer, cofactor::Integer, name::Union{Symbol, Nothing}) where P <: AbstractPoint = ECPoint{P, static_ECPoint(order, cofactor, name)}

order(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = BigInt(S.order)

cofactor(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S} = S.cofactor

name(::Type{ECPoint}) = nothing
name(::Type{ECPoint{P}}) where P <: AbstractPoint = nothing


function name(::Type{ECPoint{P, S}}) where {P <: AbstractPoint, S}
    
    str = uint2string(S.name)

    if str == ""
        return nothing
    else
        return Symbol(str)
    end
end


### May need to do epoint seperatelly
function Base.show(io::IO, ::Type{P}) where P <: ECPoint
    try
        if name(P) == nothing
            Base._show_type(io, P)
        else
            print(io, name(P))
        end
    catch e
        @warn "An odd kind of error: $e"  # UndefVarError(:P)
        print(io, "ECPoint")
    end
end

function Base.show(io::IO, p::P) where P <: ECPoint
    show(io, P)
    print(io, " <| (")
    show(io, p.p.x)
    print(io, ", ")
    show(io, p.p.y)
    print(io, ")")
end


function Base.display(::Type{P}) where P <: ECPoint
    show(P)
    ### I could be more precise on the Point
    ### Like ECPoint{AffinePoint{<:Weierstrass, <:F2GNB}, ::static_ECPoint}
    print(" (alias for ECPoint{<:AffinePoint, ::static_ECPoint})") 
end


function Base.:+(x::P, y::P) where P <: ECPoint 
    @assert x.p != y.p "AffinePoint's can't be equal. Multiply by 2 instead." 
    # Also assertion that x.p != -x.p
    return P(x.p + y.p)
end

Base.:*(x::P, n::Integer) where P <: ECPoint = P(x.p * n)
Base.:*(n::Integer, x::ECPoint) = x * n


<|(::Type{ECPoint{P, S}}, x) where {P <: AbstractPoint, S} = ECPoint{P, S}(P <| x)
validate(p::P) where P <: ECPoint = order(P) * p.p == zero(p.p) 

oncurve(p::ECPoint) = oncurve(p.p)

Base.:(==)(x::ECPoint, y::ECPoint) = x.p == y.p


gx(p::ECPoint) = gx(p.p)
gy(p::ECPoint) = gy(p.p)


