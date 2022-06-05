using Nettle
using Random: AbstractRNG
using CryptoUtils: is_quadratic_residue, sqrt_mod_prime


struct Hash
    spec::String
end

(h::Hash)(x::Vector{UInt8}) = hex2bytes(hexdigest(h.spec, x))

# Dispatching on value types seems as plausable solution
function outlen(h::Hash) 
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

#using Random: AbstractRNG

struct PRG <: AbstractRNG
    h::Hash
    s::Vector{UInt8}
end

(prg::PRG)(i::UInt32) = prg.h([prg.s..., reverse(reinterpret(UInt8, UInt32[i]))...])


function Base.getindex(prg::PRG, range)
    (; start, stop) = range
    
    a = outlen(prg.h) ÷ 8 

    K = div(stop, a, RoundUp) - 1

    r = UInt8[]
    
    for i in UInt32(0):UInt32(K)
        ri = prg(i)
        append!(r, ri)
    end
    
    return r[range]
end


struct RO
    h::Hash
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

Base.rand(prg::PRG, n::Int, N::Int) = Base.rand(prg, BigInt, N; n)



struct ROPRG
    ρ::Vector{UInt8}
    rohash::Hash
    prghash::Hash
end


function (roprg::ROPRG)(x::Vector{UInt8})

    (; ρ, rohash, prghash) = roprg

    ns = outlen(prghash)
    ro = RO(rohash, ns)

    d = UInt8[ρ..., x...]   

    s = ro(d)
    prg = PRG(prghash, s)
    return prg
end

(roprg::ROPRG)(x::String) = roprg(Vector{UInt8}(x))
(roprg::ROPRG)(x::Symbol) = roprg(string(x))







#function Base.rand(prg::PRG, ::Type{G}, N::Integer; nr::Integer = 0) where G <: PGroup

function Base.rand(prg::PRG, spec::MODP, N::Integer; nr::Integer = 0) 

    #p = modulus(G)
    #q = order(G)

    p = modulus(spec)
    q = order(spec)


    np = bitlength(p)

    𝐭 = rand(prg, BigInt, N; n = np + nr)

    𝐭′ = mod.(𝐭, big(2)^(np + nr))

    𝐡 = powermod.(𝐭′, (p - 1) ÷ q, p)
    
    #𝐡_typed = convert(Vector{PGroup{G}}, 𝐡)
    #𝐡_typed = convert(Vector{G}, 𝐡)


    #return 𝐡_typed
    return 𝐡
end



function Base.rand(prg::PRG, spec::ECP, N::Integer; nr::Integer = 0) 

    (; a, b) = spec

    p = modulus(spec)
    q = order(spec)


    np = bitlength(p) # 1

    𝐭 = rand(prg, BigInt, N*10; n = np + nr)  # OPTIMIZE

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


#<|(::Type{Vector{G}}, x::Vector) where G <: Group = G[ G <| i for i in x]


#Base.rand(prg::PRG, ::Type{G}, N::Integer; nr::Integer = 0) where G <: ECGroup = Vector{G} <| rand(prg, spec(G), N; nr)
