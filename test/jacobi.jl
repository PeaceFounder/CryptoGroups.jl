using Test
using CryptoGroups
using CryptoUtils

# Helper function to test the implementation
function test_jacobi(jacobi::Function)
    # Test cases: (x, y, expected_result)
    test_cases = [
        (1, 1, 1),
        (1, 3, 1),
        (2, 3, -1),
        (1, 5, 1),
        (2, 5, -1),
        (5, 9, 1),
        (7, 9, 1),
        (0, 3, 0),
        #(-1, 3, 1),
    ]
    
    for (x, y, expected) in test_cases
        @test jacobi(BigInt(x), BigInt(y)) == expected
    end
end

test_jacobi(CryptoUtils.jacobi)
test_jacobi(CryptoGroups.Utils.jacobi)

for i in 1:100

    x = rand(1:BigInt(2)^i)
    p = CryptoUtils.random_prime(i)

    if p == 2
        continue
    end

    @test CryptoGroups.Utils.jacobi(x, p) == CryptoUtils.jacobi(x, p)
end
