# # Lagrange Interpolation

# This Julia code demonstrates the implementation of Lagrange interpolation over a modular field, a crucial component in cryptographic schemes such as Shamir's Secret Sharing. The implementation showcases the flexibility and composability of Julia's ecosystem by seamlessly integrating the external `Polynomials` package with the custom modular field arithmetic provided by `CryptoGroups`. 

# `CryptoGroups` handles the finite field arithmetic, while `Polynomials` manages the polynomial operations, without creating a direct dependency between the two. The result is a concise yet powerful implementation that can be easily adapted for various cryptographic applications. The example includes a test case that reconstructs a secret (the constant term of the polynomial) using Lagrange interpolation, illustrating its practical application in secret sharing schemes.

using Test
using CryptoGroups.Fields
using Polynomials

function lagrange_interpolation(x::Vector{T}, y::Vector{T}) where T
    n = length(x)
    result = Polynomial{T}([0])
    
    for i in 1:n
        term = Polynomial{T}([1])
        for j in 1:n
            if i != j
                term *= Polynomial{T}([-x[j], 1]) / (x[i] - x[j])
            end
        end
        result += y[i] * term
    end
    
    return result
end

p = 23 
secret = 3

poly = Polynomial{FP{p}}([secret, 5, 2, 5])

x = FP{p}[2, 4, 6, 8]
y = poly.(x)

interp_poly = lagrange_interpolation(x, y)

@test value(interp_poly(0)) == secret
