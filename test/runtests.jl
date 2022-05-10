using Test

@testset "Galois field tests" begin
    include("galois_fields.jl")
end

@testset "Elliptic curve primitive tests" begin
    include("elliptic_curves.jl")
end

@testset "Testing NIST elliptic curves" begin
    include("curve_specs.jl")
end

@testset "Testing external fields" begin
    include("../examples/external_fields.jl")
end

@testset "Testing groups and ElGamal" begin
    include("groups.jl")
    include("elgamal.jl")
end

@testset "Testing independent basis generation" begin
    include("primitives.jl")
    include("gbasis.jl")
    include("gecbasis.jl")
end


