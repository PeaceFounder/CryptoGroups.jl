# # Lagrange Interpolation
# Lagrange interpolation for modular field. In practice it is used for 
# Shamir Lagrange secret sharing scheme where the constant term is the secret

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
