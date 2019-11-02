using Mods
using Primes

struct PrimeGroup <: AbstractGroup 
    G::Mod
    t::Int # security of the group
    function PrimeGroup(G,t)
        @assert isprime(G.mod)
        # I could also test by calculating Euler totient
        new(G,t)
    end
end

PrimeGroup(g,p,t) = PrimeGroup(Mod(g,p),t)
value(G::PrimeGroup) = G.G.val
security(G::PrimeGroup) = G.t


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
            return PrimeGroup(g,p,t)
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
            return PrimeGroup(g,tp)
        end
    end
end
