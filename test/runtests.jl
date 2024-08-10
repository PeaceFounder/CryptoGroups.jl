using Test

@testset "Galois field tests" begin
    include("galois_fields.jl")
end

@testset "Elliptic curve primitive tests" begin
    include("elliptic_curves.jl")
end

@testset "Testing NIST elliptic curves" begin
    #include("curve_specs.jl")
end

@testset "Testing external fields" begin
    include("../examples/external_fields.jl")
end

@testset "Testing CSPRG" begin
    include("primitives.jl")
end

@testset "Testing CryptoProofs" begin
    include("gbasis.jl")
    include("gecbasis.jl")
    include("elgamal.jl")
end

@testset "Testing Specs module" begin
    include("spec.jl")
    include("ecdsa/ecdsa_prime.jl")
    include("ecdsa/ecdsa_pb.jl")
end

@testset "Testing Conversions" begin
    include("conversions.jl")
end

@testset "Testing group API" begin
    include("groups.jl")
end



