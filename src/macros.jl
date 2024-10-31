using .Specs: modp_spec

macro PGroup(expr)
    if expr.head == :braces
        if length(expr.args) == 1 && !(expr.args[1] isa Expr)
            # Single argument case: @PGroup{some_name}
            name = expr.args[1]
            spec = modp_spec(name)
            group = concretize_type(PGroup, spec; name)
            return group

        else
            # Two-argument case: @PGroup{p=23, q=11} or @PGroup{p=my_p, q=my_q}
            p = q = nothing
            for arg in expr.args
                if arg isa Expr && arg.head == :(=)
                    lhs = arg.args[1]
                    rhs = arg.args[2]
                    if lhs == :p
                        p = rhs
                    elseif lhs == :q
                        q = rhs
                    end
                end
            end

            # Check if both p and q are provided
            if isnothing(p) || isnothing(q)
                error("Both p and q must be specified in @PGroup{p=..., q=...}")
            end

            # Properly escape both p and q values
            return quote
                local p_val = $(esc(p))
                local q_val = $(esc(q))
                local spec = MODP(p_val; q=q_val)
                concretize_type(
                    PGroup,
                    spec
                )
            end
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
    if expr.head == :braces && length(expr.args) == 1
        arg = expr.args[1]
        # Create a GlobalRef to @ECPoint in the current module
        #ecpoint_macro = GlobalRef(__module__, Symbol("@ECPoint"))

        ecpoint_macro = getproperty(@__MODULE__, Symbol("@ECPoint"))

        point_expr = Expr(:macrocall, ecpoint_macro, LineNumberNode(@__LINE__), 
                         Expr(:braces, arg))

        return :($ECGroup{$(esc(point_expr))})
    else
        error("Invalid syntax. Use @ECGroup{curve_name} or @ECGroup{Module.curve_name}")
    end
end

# First, let's modify @ECPoint to ensure it handles symbol quoting correctly
macro ECPoint(expr)
    if expr.head == :braces && length(expr.args) == 1
        arg = expr.args[1]
        
        # Handle module-qualified names (e.g., OpenSSLGroups.SecP256k1)
        if arg isa Expr && arg.head == :. 
            return quote
                local P = $(esc(arg))
                concretize_type(ECPoint{P}, order(P), cofactor(P); name = nameof(P))
            end
        # Handle simple symbols
        elseif arg isa Symbol
            # Important: Use QuoteNode here for the isdefined check
            return quote
                if isdefined($(__module__), $(QuoteNode(arg)))
                    # If defined, use the escaped symbol to access its value
                    local P = $(esc(arg))
                    concretize_type(ECPoint{P}, order(P), cofactor(P); name = nameof(P))
                else
                    # If not defined, treat it as a curve name
                    local spec = curve($(QuoteNode(arg)))
                    concretize_type(ECPoint, spec)
                end
            end
        else
            error("Invalid syntax. Use @ECPoint{curve_name} or @ECPoint{Module.curve_name}")
        end
    else
        error("Invalid syntax. Use @ECPoint{curve_name} or @ECPoint{Module.curve_name}")
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
