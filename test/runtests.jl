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

@safetestset "Testing Specs module" begin
    include("spec.jl")
    include("ecdsa/ecdsa_prime.jl")
    include("ecdsa/ecdsa_pb.jl")
end

@safetestset "Testing Conversions" begin
    include("conversions.jl")
end

@safetestset "Testing Jacoby symbol" begin
    include("jacobi.jl")
end

@safetestset "Testing group API" begin
    include("groups.jl")
end

### Examples
println("\n\nProceeding with examples\n")

@safetestset "External Fields" begin
    include("../examples/external_fields.jl")
end

@safetestset "Reed Solomon Error Correction" begin
    include("../examples/reed-solomon.jl")
end

@safetestset "Lagrange Over Prime Field" begin
    include("../examples/lagrange.jl")
end

@safetestset "Digital Signature Algorithm" begin
    include("../examples/dsa.jl")
end

@safetestset "Key Encapsulation Mechanism" begin
    include("../examples/kem.jl")
end

@safetestset "ElGamal Cryptosystem" begin
    include("../examples/elgamal.jl")
end

@safetestset "Proof of Knowledge" begin
    include("../examples/knowledge_proof.jl")
end


