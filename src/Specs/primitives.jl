using Nettle
import Random
using Random: AbstractRNG

using CryptoUtils: is_quadratic_residue, sqrt_mod_prime

struct HashSpec
    spec::String
end

(h::HashSpec)(x::Vector{UInt8}) = hex2bytes(hexdigest(h.spec, x))
(h::HashSpec)(x::String) = h(Vector{UInt8}(x))


# Dispatching on value types seems as plausable solution
function bitlength(h::HashSpec) 
    s = h.spec

    if s == "sha256"
        return 256
    elseif s == "sha384"
        return 384
    elseif s == "sha512"
        return 512
    else
        error("No corepsonding mapping for $x implemented")
    end
end

struct PRG <: AbstractRNG
    h::HashSpec
    s::Vector{UInt8}
end

PRG(hasher::String; s = Vector{UInt8}("SEED")) = PRG(HashSpec(hasher), s)
bitlength(prg::PRG) = bitlength(prg.h)

(prg::PRG)(i::UInt32) = prg.h([prg.s..., reverse(reinterpret(UInt8, UInt32[i]))...])


function Base.getindex(prg::PRG, range)
    (; start, stop) = range
    
    a = bitlength(prg.h) ÷ 8 # outlength

    K = div(stop, a, RoundUp) - 1

    r = UInt8[]
    
    for i in UInt32(0):UInt32(K)
        ri = prg(i)
        append!(r, ri)
    end
    
    return r[range]
end


struct RO
    h::HashSpec
    n_out::Int
end

zerofirst(x, n) = (x << n) >> n # Puts first n bits of a number x to zero. 

function (ro::RO)(d::Vector{UInt8})
    (; h, n_out) = ro

    nb = reinterpret(UInt8, UInt32[n_out])
    s = h([reverse(nb)...,d...]) # Numbers on Java are represented in reverse
    prg = PRG(h, s)

    a = prg[1:div(n_out, 8, RoundUp)]
    
    if mod(n_out, 8) != 0 
        a[1] = zerofirst(a[1], 8 - mod(n_out, 8))
    end

    return a
end

_tobig(x) = parse(BigInt, bytes2hex(reverse(x)), base=16)
interpret(::Type{BigInt}, x::Vector{UInt8}) = _tobig(reverse(x))


function interpret(::Type{Vector{T}}, 𝐫::Vector{UInt8}, N::Int) where T <: Integer
    M = length(𝐫) ÷ N
    𝐮 = reshape(𝐫, (M, N))
    𝐭 = [interpret(T, 𝐮[:, i]) for i in 1:N]
    return 𝐭
end



function Base.rand(prg::PRG, ::Type{T}, N::Int; n = bitlength(T)) where T <: Integer

    M = div(n, 8, RoundUp) # bytes for each number

    total = M * N

    𝐫 = prg[1:total]
    𝐭 = interpret(Vector{BigInt}, 𝐫, N)
    
    return 𝐭
end

#Base.rand(prg::PRG, n::Int, N::Int) = Base.rand(prg, BigInt, N; n)
@deprecate Base.rand(prg::PRG, n::Int, N::Int) Base.rand(prg, BigInt, N; n) false
Base.rand(prg::PRG, ::Type{T}; n = bitlength(T)) where T <: Integer = rand(prg, T, 1; n)[1]

function Random.rand!(rng::PRG, a::AbstractArray{T}, sp) where T <: Integer

    values = rand(rng, BigInt, length(a); n = bitlength(maximum(sp))) 

    a_flat = reshape(a, length(a))
    a_flat .= minimum(sp) .+ mod.(values, maximum(sp) - minimum(sp) + 1) # ToDo: fix bias (a simple skipping strategy should work)

    return a
end


struct ROPRG
    ρ::Vector{UInt8}
    rohash::HashSpec
    prghash::HashSpec
end


function (roprg::ROPRG)(x::Vector{UInt8})

    (; ρ, rohash, prghash) = roprg

    ns = bitlength(prghash) # outlen
    ro = RO(rohash, ns)

    d = UInt8[ρ..., x...]   

    s = ro(d)
    prg = PRG(prghash, s)
    return prg
end

(roprg::ROPRG)(x::String) = roprg(Vector{UInt8}(x))
(roprg::ROPRG)(x::Symbol) = roprg(string(x))


function modp_generator_basis(prg::PRG, p::Integer, q::Integer, N::Integer; nr::Integer = 0)

    np = bitlength(p)

    𝐭 = rand(prg, BigInt, N; n = np + nr)

    𝐭′ = mod.(𝐭, big(2)^(np + nr))

    𝐡 = powermod.(𝐭′, (p - 1) ÷ q, p)
    
    return 𝐡
end

function ecp_generator_basis(prg::PRG, (a, b)::Tuple{Integer, Integer}, p::Integer, q::Integer, N::Integer; nr::Integer = 0)

    np = bitlength(p) # 1

    𝐭 = rand(prg, BigInt, N*10; n = np + nr)  # OPTIMIZE (I would need it as an iterator)

    𝐭′ = mod.(𝐭, big(2)^(np + nr))

    𝐳 = mod.(𝐭′, p)

    𝐡 = Vector{Tuple{BigInt, BigInt}}(undef, N)

    l = 1

    f(x) = x^3 + a*x + b # This assumes that I do know how to do arithmetics with fields.

    for zi in 𝐳
        y2 = mod(f(zi), p)

        if is_quadratic_residue(y2, p)

            x = zi
            y = sqrt_mod_prime(y2, p)

            # The smallest root is taken
            if p - y < y
                y = p - y
            end

            𝐡[l] = (x, y)

            if l == N
                break
            else
                l += 1                
            end
        end
    end

    if l != N
        error("Not enough numbers for 𝐭 have been allocated")
    end

    return 𝐡
end

# ToDo: consider deprecating and redirect to generator_basis function
function Base.rand(prg::PRG, spec::MODP, N::Integer; nr::Integer = 0) 

    p = modulus(spec)
    q = order(spec)

    @assert !isnothing(q) "Order of the group must be known"

    return modp_generator_basis(prg, p, q, N; nr)
end

Base.rand(prg::PRG, spec::ECP, N::Integer; nr::Integer = 0) = ecp_generator_basis(prg, (spec.a, spec.b), modulus(spec), order(spec), N; nr)

