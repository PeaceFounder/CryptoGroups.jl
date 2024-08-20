# CryptoGroups

Cryptographic groups are a fundamental building block for digital signatures, key exhange algorithm, assymetric encryption and many other exciting algorithms of practical importance. 

![](https://raw.githubusercontent.com/PeaceFounder/CryptoGroups.jl/b7e6d4b8be1807e124422229428bb4c289523769/doc/assets/CryptoGroups%20types.svg) 

## Elliptic Curve safety

- Elliptic curve safety is implemented in `ECPoint` wrapper type that adds essential safety checks;
- Elliptic curve point is rejected at `ECPoint` constructor if `p^h == 0` or if it does not satisfy the corresponding curve equation (see `oncurve` method);
- According to https://safecurves.cr.yp.to/complete.html, "The standard Weierstrass addition formulas fail if Q happens to match -P" is checked; 
