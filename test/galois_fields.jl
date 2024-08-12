using Test
import CryptoGroups.Fields: F2PB, F2GNB, print_poly, red!, mul, mul_gnb, construct_integer_order_prime, FP, order, reducer, tobits
import CryptoGroups.Utils: @bin_str


f = bin"10011"
c = reverse(bin"1100101")
@test red!(c, reverse(f)) == bin"1111" # I already changed rem and order functions and that this is the result

### Additional test for multiplication

@test mul(bin"0100", bin"0100") == bin"00100000"

a = reverse(bin"1101")
b = reverse(bin"1001")

c = mul(a, b)

@test red!(c, reverse(f)) == bin"1111"

@test mul(reverse(bin"1101"), reverse(bin"1001"), reverse(bin"10011")) == bin"1111"
@test mul(reverse(bin"0010"), reverse(bin"0010"), reverse(bin"10011")) == reverse(bin"0100")
@test mul(reverse(bin"0010"), reverse(bin"0100"), reverse(bin"10011")) == reverse(bin"1000")

a = F2PB(f)(bin"1101")
b = F2PB(f)(bin"1001")

c = F2PB(f)(bin"1111")

@test a*b == c

### Testing generation of basis elements

α = F2PB(f)(bin"0010")

@test α^(order(α) + 1) == α

@test α^8 * α^8 == α^16
@test α * α * α * α * α == α^5
@test α^14 * α == α^15

@test α^3 * α^5 == α^8

@test inv(α) * α == one(α)


@test zero(F2GNB{5, 1}) == F2GNB{5, 1}(bin"00000")
@test one(F2GNB{5, 1}) == F2GNB{5, 1}(bin"11111")

a = F2GNB{4, 1}(bin"1101")
b = F2GNB{4, 1}(bin"1011")

@test a + b == F2GNB{4, 1}(bin"0110")


@test construct_integer_order_prime(3, 3*4 + 1) in [3, 9, 1]

### Now the last bit of getting multiplication right 

@test mul_gnb(bin"1000", bin"1101", 3) == bin"0010"

### 

a = F2GNB{4, 3}(bin"1000")
b = F2GNB{4, 3}(bin"1101")

c = F2GNB{4, 3}(bin"0010")

@test a * b == c

@test c^(order(c) + 1) == c
@test inv(c) * c == one(c)

a = FP{23}(17)
b = FP{23}(3)

@test a * b == FP{23}(5)

@test b^11 == FP{23}(1)


#### More diligent testing for polynomial basis

# function field2poly(z::BitVector)
#     poly = Int[]
#     for (i, b) in enumerate(z)
#         if b==true
#             push!(poly, i - 1)
#         end
#     end
#     return poly
# end

#field2poly(z) = field2poly(z.x)

# Arguments are written with little endianess in mind
# let


#     F = F2PB([163, 7, 6, 3, 0])
#     f = reducer(F)

#     x19 = F(BitVector(i==20 for i in 1:163))
#     @test findfirst(x -> x== 1, mul(x19.x, x19.x)) == 39

#     x38 = F(BitVector(i==39 for i in 1:163))
#     @test findfirst(x -> x==1, red!(x38.x, f)) == 39

#     @test field2poly(x19 * x19) == [38]

#     x144 = F(BitVector(i==145 for i in 1:163))
#     x145 = F(BitVector(i==146 for i in 1:163))
#     @test field2poly(x145 * x144) == [126, 129, 132, 133]

#     x160 = F(BitVector(i==161 for i in 1:163))
#     x159 = F(BitVector(i==160 for i in 1:163))

#     @test field2poly(x160 * x159) == [0, 3, 6, 7, 156, 159, 162] 

#     _x319 = BitVector(i==320 for i in 1:163*2)
#     @test field2poly(red!(_x319, f)) == [0, 3, 6, 7, 156, 159, 162]

# end
