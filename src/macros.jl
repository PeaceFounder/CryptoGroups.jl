using .Specs: modp_spec


# TODO: Add support for @PGroup{p = _p, q = _q} where _q and _p are defined out of the scope
macro PGroup(expr)
    if expr.head == :braces
        if length(expr.args) == 1 && !(expr.args[1] isa Expr)
            # Single argument case: @PGroup{some_name}
            name = expr.args[1]
            spec = modp_spec(name)
            group = concretize_type(PGroup, spec; name)
            return group
        else
            # Two-argument case: @PGroup{p=23, q=11}
            p = q = nothing
            for arg in expr.args
                if arg isa Expr && arg.head == :(=)
                    if arg.args[1] == :p
                        p = arg.args[2]
                    elseif arg.args[1] == :q
                        q = arg.args[2]
                    end
                end
            end

            # Check if both p and q are provided
            if isnothing(p) || isnothing(q)
                error("Both p and q must be specified in @PGroup{p=..., q=...}")
            end

            spec = MODP(p; q)
            group = concretize_type(PGroup, spec)
            return group
        end
    else
        error("Invalid syntax. Use @PGroup{p=..., q=...} or @PGroup{some_name}")
    end
end


Base.show(io::IO, ::Type{PGroup}) = print(io, "PGroup")

function Base.show(io::IO, ::Type{G}) where G <: PGroup
    if @isdefined G

        print(io, "@PGroup{")

        if name(G) != nothing

            print(io, string(name(G)))          

        else
            
            p = modulus(G)
            q = order(G)

            print(io, "p = $p, q = $q")
            
        end

        print(io, "}")

    else
        print(io, "PGroup")
    end
end


macro ECGroup(expr)
    if expr.head == :braces && length(expr.args) == 1 && !(expr.args[1] isa Expr)
        # Single argument case: @PGroup{some_name}
        some_name = expr.args[1] 

        # If the curve can't be found error here
        spec = curve(some_name)
        group = concretize_type(ECGroup, spec)

        return group
    else
        error("Invalid syntax. Use @ECGroup{curve_name}")
    end
end

function Base.show(io::IO, g::G) where G <: ECGroup
    show(io, G)
    print(io, "(")
    show(io, gx(g))
    print(io, ", ")
    show(io, gy(g))
    print(io, ")")
end

function Base.show(io::IO, ::Type{G}) where G <: ECGroup
    if @isdefined G
        if name(G) == nothing
            Base._show_type(io, G)
        else
            print(io, "@ECGroup{$(name(G))}")
        end
    else
        print(io, "ECGroup")
    end
end


macro ECPoint(expr)
    if expr.head == :braces && length(expr.args) == 1 && !(expr.args[1] isa Expr)
        # Single argument case: @PGroup{some_name}
        some_name = expr.args[1] 

        spec = curve(some_name)
        # If the curve can't be found error here
        _curve = concretize_type(ECPoint, spec)

        return _curve
    else
        error("Invalid syntax. Use @ECPoint{curve_name}")
    end
end


### May need to do epoint seperatelly
function Base.show(io::IO, ::Type{P}) where P <: ECPoint
    if @isdefined P
        if name(P) == nothing
            Base._show_type(io, P)
        else
            print(io, "@ECPoint{$(name(P))}")
        end
    else
        print(io, "ECPoint")
    end
end


function Base.show(io::IO, p::P) where P <: ECPoint
    show(io, P)
    print(io, "(")
    show(io, gx(p))
    print(io, ", ")
    show(io, gy(p))
    print(io, ")")
end

function Base.display(::Type{P}) where P <: ECPoint
    show(P)
end



