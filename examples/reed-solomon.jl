# # Reed-Solomon Error Correction
# This is Reed-Solomon algorithm implemnetd in Julia. It demonstraates that polynomial field
# ease of use and ability to compose it with Polynomials library which can do polynom multiplication, 
# derivatives and modular reduction.

# HELP WANTED: Currently the implementation has bugs and can't even get correct error positions; To proceed
# one could write a tests for `berklekamp_massey` and when it works debug error_positions

using Test
using CryptoGroups.Fields
using Polynomials

GF256 = @F2PB{X^8 + X^4 + X^3 + X + 1}

roots(p) = GF256[i for i in 0:255 if iszero(p(GF256(i)))]

function berklekamp_massey(syndromes::Vector{GF256})
    B = Polynomial{GF256}([1])
    C = Polynomial{GF256}([1])
    
    L, m = 0, 1
    for n in 1:length(syndromes)
        ## Adjust for zero-based indexing in Polynomial
        d = syndromes[n] + sum(C[i-1] * syndromes[n - i + 1] for i in 1:L; init = zero(GF256))
        if iszero(d)
            m += 1
        elseif 2*L <= n
            T = C
            ## Use zero-based indexing for polynomial multiplication
            C -= d * B * Polynomial{GF256}([0, 1])^(m-1)
            L = n + 1 - L
            B = T
            m = 1
        else
            ## Use zero-based indexing for polynomial multiplication
            C -= d * B * Polynomial{GF256}([0, 1])^(m-1)
            m += 1
        end
    end
    
    return C
end

function rs_encode(message::Vector{UInt8}, g::Polynomial{GF256})

    m = Polynomial{GF256}(message)

    r_poly = m % g
    r_octets = r_poly[:] .|> (x -> octet(x)[1])

    return [message..., r_octets...]
end

function rs_decode(received::Vector{UInt8}, n::Int, k::Int, g::Polynomial{GF256})
    
    received_poly = Polynomial{GF256}(received)

    syndromes = [received_poly(r) for r in roots(g)]

    if all(iszero(s) for s in syndromes)
        return received[1:k] ## no errors detected
    end
    
    error_locator = berklekamp_massey(syndromes)

    error_positions = [i for i in 1:n if iszero(error_locator(GF256(2)^i))]
    ## @show error_positions

    error_evaluator = Polynomial{GF256}(syndromes) * error_locator % Polynomial{GF256}([zeros(Int, n - k)..., 1])

    errors = zeros(GF256, n)

    for i in error_positions
        x = GF256(2) ^ (n - i)
        y = error_evaluator(x) / derivative(error_locator)(x)
        errors[i] = y
    end

    corrected = GF256[received...] .- errors

    return corrected[1:k] .|> (x -> octet(x)[1])
end

n, k = 15, 11
g = Polynomial{GF256}([1, 25, 6, 8, 14])

message = rand(UInt8, k)
encoded = rs_encode(message, g)

errors = zeros(UInt8, n)
errors[3] = rand(0:255)
errors[7] = rand(0:255)
received = encoded + errors

decoded = rs_decode(received, n, k, g)

## @test message == decoded
