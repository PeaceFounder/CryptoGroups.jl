using Mods
using Primes

struct PrimeGroup <: CyclicGroup 
    G::Mod
    q::Integer
    t::Int # security of the group
    function PrimeGroup(G,q,t)
        @assert isprime(G.mod)
        # I could also test by calculating Euler totient
        new(G,q,t)
    end
end

PrimeGroup(g,p,q,t) = PrimeGroup(Mod(g,p),q,t)
PrimeGroup(val::Integer,G::PrimeGroup) = PrimeGroup(val,G.G.mod,G.q,G.t)

value(G::PrimeGroup) = G.G.val
order(G::PrimeGroup) = G.q
security(G::PrimeGroup) = G.t

import Base.*
function *(X::PrimeGroup,Y::PrimeGroup)
    # Parametrizing group with p and q is not an option as that would recompile the code. 
    qx,qy = order(X), order(Y)
    px,py = X.G.mod, Y.G.mod
    
    if px!=py || qx!=qy
        error("Groups are not equal")
    else
        G = X.G * Y.G
        return PrimeGroup(G,X.q,X.t)
    end
end

import Base.==
function ==(X::PrimeGroup,Y::PrimeGroup)
    qx,qy = order(X), order(Y)
    px,py = X.G.mod, Y.G.mod
    
    if px!=py || qx!=qy
        error("Groups are not equal")
    else
        return X.G==Y.G
    end
end

### Some prime group generation algorithms. References:
# + https://crypto.stackexchange.com/questions/820/how-does-one-calculate-a-primitive-root-for-diffie-hellman
# + https://math.stackexchange.com/questions/124408/finding-a-primitive-root-of-a-prime-number

### I could manage to get 512 security
function SophieGermainGroup(rng::AbstractRNG,g::Integer,t::Int)
    @assert g!=1
    while true
        q = rngprime(rng,2*t)
        p = 2*q + 1
        if g!=p-1 && isprime(p)
            return PrimeGroup(g,p,2*q,t)
        end
    end
end

function generateqp(rng::AbstractRNG,t)
    while true
        q = rngprime(rng,t)
        r = rngint(rng,t)

        p = q*r + 1
        if isprime(p)
            return q,p
        end
    end
end

function DSAStandartGroup(rng::AbstractRNG,tp::Int,tg::Int)
    q,p = generateqp(rng,2*tp)
    while true
        u = Mod(rngint(rng,tg),p)
        g = u^div((p-1),q)
        if g.val!=1 || g.val!=0
            return PrimeGroup(g,q,tp)
        end
    end
end
