# CryptoGroups
[![Build Status](https://travis-ci.com/PeaceFounder/CryptoGroups.jl.svg?branch=master)](https://travis-ci.com/PeaceFounder/CryptoGroups.jl)

Cryptographic groups are a fundamental building block for digital signatures, key exhange algorithm, assymetric encryption and many other exciting algorithms of practical importance. 

![](https://raw.githubusercontent.com/PeaceFounder/CryptoGroups.jl/b7e6d4b8be1807e124422229428bb4c289523769/doc/assets/CryptoGroups%20types.svg) 

## ToDo

  * [x] Import and fix tests
  * [x] Add a `spec` function with which specs can be retrieved as `spec(:P_192)`, `spec(:OakleyV1)` or `spec(:B_163, :PB)`.
  * [x] Introuce abstract type `Spec`
  * [x] Rename `crs` to `rand` and in `ShuffleProofs`, `gen_verificatum_prg`.
  * [x] Rename `solidify` as `specialize`
  * [x] Rename `incurve` to `oncurve`
  * [ ] According to https://safecurves.cr.yp.to/complete.html 

    > ... the standard Weierstrass addition formulas fail if Q happens to match -P. This will not be caught by random tests. 
  
  as well as identical points can not be summed. Could be partially addressed at the higher level of `ECGroup`.

    > An implementor can stop a small-subgroup attack by rejecting any Q for which hQ = 0
  
  This may be addressed at constructor level, but requires to know the cofactor. 

  * [x] Adding accessor methods to `AffinePoint` as `_a` and `_b` and acessor methods to curves `a` and `b`
  * [x] Implement independent basis generation for elliptic curves
      * [x] Add a square root function for elliptic curves (Imported from CryptoUtils)
  * [ ] Make a prg iterator for numbers
  * [x] Fix the `UndefVarError(:P)` in the show method
  * [x] Add point, field, integer conversions as specified in X9.62 section 4.3 
  * [ ] Specify cofactors in the elliptic curve specs and encode cofactor assertions in `ECPoint`
  * [ ] Does order needs to be computed from `n` by divifing with cofactor `h`?
  * [ ] Add some docs
  * [ ] Consider better alternatives for internal data representation of `F2GNB` and `F2PB` to improve performance.
