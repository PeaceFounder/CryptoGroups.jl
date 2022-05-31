struct FP{P} <: PrimeField
    x::BigInt # Instead of BigInt I could use BitIntegers here. Need to check performance of powermod...
    FP{P}(x::Integer) where P = new{P}(BigInt(x))
end

specialize(::Type{FP}, p::Integer) = FP{StaticBigInt(p)} # A workaround to store BigInt as type argument

value(x::FP) = x.x
modulus(::FP{P}) where P = BigInt(P)
modulus(::Type{FP{P}}) where P = BigInt(P)
order(::FP{P}) where P = BigInt(P)
order(::Type{FP{P}}) where P = BigInt(P)


function trimnumber(x::String)
    if length(x) < 30
        return x
    else
        return x[1:10] * "..." * x[end-10:end]
    end
end

trimnumber(x::Integer)= trimnumber(string(x))

modulus_str(x::BigInt) = trimnumber(x)
modulus_str(::Type{F}) where F <: FP = modulus_str(modulus(F))

Base.show(io::IO, ::Type{FP}) = print(io, "FP")
Base.show(io::IO, ::Type{F}) where F <: FP = print(io, "FP/$(modulus_str(F))")

#########################################

### StaticBitVEctor would be a better fit here

struct Reducer{N}
    len::Int
    chunks::NTuple{N, UInt64}
end


Reducer(x::BitVector) = Reducer(x.len, Tuple(x.chunks))

Reducer(x::Vector{Int}) = Reducer(BitVector(i in x for i in 0:maximum(x)))

Base.length(x::Reducer) = x.len

binary_order(R::Reducer) = length(R) - 1


function Base.convert(::Type{BitVector}, x::Reducer)

    a = BitVector(undef, x.len)
    a.chunks = UInt64[x.chunks...]

    return a
end

function xstr(i::Int)
    if i == 0
        return "1"
    elseif i == 1
        return "X"
    else
        return "X^$i"
    end
end

get_powers(x::BitVector) = (0:x.len - 1)[x]
get_powers(x::Reducer) = get_powers(convert(BitVector, x))

function print_poly(io::IO, x::BitVector)
    n = reverse(get_powers(x)) # Best to visualize dominant first
    write(io, "$(xstr(n[1]))")
    for i in n[2:end]
        write(io, " + $(xstr(i))")
    end
end

print_poly(io::IO, x::Reducer) = print_poly(io, convert(BitVector, x)) 
print_poly(x::Any) = print_poly(Base.stdout, x)

Base.show(io::IO, x::Reducer) = print_poly(io, x)


struct F2PB{R} <: BinaryField ### R is reducer
    x::BitVector

    function F2PB{R}(x::BitVector) where R
        @assert length(x) == length(R) - 1
        new(x)
    end
end


F2PB(x::F2PB) = x

specialize(::Type{F2PB}, f::Vector{Int}) = F2PB{Reducer(f)}


binary_order(::Type{F2PB{R}}) where R = length(R) - 1
binary_order(::F) where F <: F2PB = binary_order(F)




function frombits(::Type{F}, bits::BitVector) where F <: F2PB
    m = binary_order(F)
    truncated = reverse(bits)[1:m] 
    return F(truncated)
end

tobits(x::F2PB) = reverse(x.x)


Base.zero(::Type{F2PB{R}}) where R = F2PB{R}(BitVector(0 for i in 1:length(R)-1))
Base.one(::Type{F2PB{R}}) where R = F2PB{R}(BitVector((1, (0 for i in 2:length(R)-1 )...)))

Base.:+(x::F, y::F) where F <: F2PB = F(xor.(x.x, y.x))

order(::Type{F}) where F <: F2PB = BigInt(2)^binary_order(F) - 1
order(x::F) where F <: F2PB = order(F)    


_order(x::BitVector) = findlast(x->x==1, x) - 1
_rem(x::BitVector) = x[1:findlast(x->x==1, x)-1]

function red!(c::BitVector, f::BitVector)

    m = _order(f)

    r = _rem(f)
    append!(r, BitVector(0 for i in 1:m))

    s = 2m - length(c)
    
    append!(c, BitVector(0 for i in 1:s))
    
    for i in length(c):-1:m+1 ### iteration may be reverse order
        if c[i] == true
            ri = r >> ((i - 1) - m )
            c = xor.(c, ri)
        end
    end

    return c[1:m]
end


function mul(a::BitVector, b::BitVector)
    
    a = copy(a)

    n = length(a)

    append!(a, BitVector(0 for i in 1:n))

    s = zero(a)

    for i in 1:n
        if b[i]
            s .⊻= a >> (i - 1)
        end
    end

    return s
end

mul(a::BitVector, b::BitVector, f::BitVector) = red!(mul(a, b), f)
Base.:*(x::F2PB{R}, y::F2PB{R}) where R = F2PB{R}(mul(x.x, y.x, convert(BitVector, R)))
Base.:(==)(x::F, y::F) where F <: F2PB = x.x == y.x

#####################################################################################

struct F2GNB{N, T} <: BinaryField

    x::BitVector

    function F2GNB{N, T}(x::BitVector) where {N, T}
        @assert length(x) == N
        new(x)
    end

end

specialize(::Type{F2GNB}, N::Int, T::Int) = F2GNB{N, T} 

### Theese two methods are essential for the field so that a field element would be possible to be iitialized 
### from external sources
frombits(::Type{F2GNB{N, T}}, bits::BitVector) where {N, T} = F2GNB{N, T}(bits[end - N + 1:end])
tobits(x::F2GNB) = x.x

binary_order(::Type{F2GNB{N, T}}) where {N, T} = N
binary_order(::F) where F <: F2GNB = binary_order(F)

order(::Type{F}) where F <: F2GNB = BigInt(2)^binary_order(F) - 1
order(x::F) where F <: F2GNB = order(F)    


Base.zero(::Type{F2GNB{N, T}}) where {N, T} = F2GNB{N, T}(BitVector(false for i in 1:N))
Base.one(::Type{F2GNB{N, T}}) where {N, T} = F2GNB{N, T}(BitVector(true for i in 1:N))


square(x::F) where F <: F2GNB = F(circshift(x.x, 1))


Base.:+(x::F, y::F) where F <: F2GNB = F(xor.(x.x, y.x))
Base.:(==)(x::F, y::F) where F <: F2GNB = x.x == y.x


# Improving performace of this function is top priority
# Random access order with Fvec is likelly making a lot of chache misses

function first_coordinate_product(a::BitVector, b::BitVector, T::P, F::Vector{P}) where P <: Integer

    m = length(a)

    p = T*m + 1

    J = 0
    
    if mod(T, 2) != 0

        # In the spec it is given that summation index goes to m/2
        # If m/2 would be performed as ceiling 
        # The following code would ask for an element like b[m + 1]
        # Alternativelly a leading bit could be added then to make araray even
        for k in 1:div(m, 2)  
            J ⊻= a[k] * b[div(m, 2) + k]
            J ⊻= a[div(m, 2) + k] * b[k]
        end
    end

    c0 = J

    for k in 1:(p - 2)
        c0 ⊻= a[F[k + 1] + 1] * b[F[p - k] + 1]
    end
    
    return c0

end

# Is a dublicate, I could import it from Specs when dust settles
function compute_integer_order(g::T, p::T) where T <: Integer
    
    @assert 1 < g < p

    b = g

    j = 1 ### Could there be an error?

    while b > 1
        b = mod(g*b, p)
        j += 1
    end

    return j
end

function construct_integer_order_prime(T, p)

    for g in 2:(p-1)

        k = compute_integer_order(g, p)

        if mod(k, T) == 0 
            u = powermod(g, div(k, T), p)
            return u
        end
    end

    error("No success")
end


function compute_F(m::Int, T::Int)
    
    p = T*m + 1

    u = construct_integer_order_prime(T, p)

    F = Vector{Int}(undef, p)

    w = 1

    for j in 0:(T - 1)
        
        n = w
        
        for i in 0:(m - 1)
            F[n] = i
            n = mod(2n, p)
        end

        w = mod(u*w, p)
    end
    

    return F
end

function mul_gnb(a::BitVector, b::BitVector, T::Int)

    u = copy(a)
    v = copy(b)

    m = length(a)

    F = compute_F(m, T)

    c = zero(u)

    # botleneck here
    for k in 0:(m - 1)
        c[k+1] = first_coordinate_product(u, v, T, F)
        circshift!(u, -1)
        circshift!(v, -1)
    end

    return c
end

Base.:*(x::F2GNB{N, T}, y::F2GNB{N, T}) where {N, T} = F2GNB{N, T}(mul_gnb(x.x, y.x, T)) 
