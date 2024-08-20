using SafeTestsets

@safetestset "Galois field tests" begin
     include("galois_fields.jl")
end

@safetestset "Elliptic curve primitive tests" begin
    include("elliptic_curves.jl")
end

@safetestset "Testing NIST elliptic curves" begin
    include("curve_specs.jl")
end

@safetestset "Testing external fields" begin
    include("../examples/external_fields.jl")
end

@safetestset "Testing Specs module" begin
    include("spec.jl")
    include("ecdsa/ecdsa_prime.jl")
    include("ecdsa/ecdsa_pb.jl")
end

@safetestset "Testing Conversions" begin
    include("conversions.jl")
end

@safetestset "Testing group API" begin
    include("groups.jl")
end



