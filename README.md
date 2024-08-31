# CryptoGroups

[![codecov](https://codecov.io/gh/PeaceFounder/CryptoGroups.jl/graph/badge.svg?token=G9HT9VSV4T)](https://codecov.io/gh/PeaceFounder/CryptoGroups.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PeaceFounder.github.io/CryptoGroups.jl/dev)

Cryptographic groups are a fundamental building block for digital signatures, key exchange algorithms, asymmetric encryption and many other exciting algorithms of practical importance. 

![](https://raw.githubusercontent.com/PeaceFounder/CryptoGroups.jl/b7e6d4b8be1807e124422229428bb4c289523769/docs/src/assets/CryptoGroups%20types.svg) 

- Cryptographic group library with generic elliptic curve implementation. It avoids the need for partial initialisations and uses the Julia type system with extensive use of type parameters, which are used for group-specifying constants. This offers unprecedented flexibility for implementing new curves and ensures the highest grade type safety. The typing is more static than what is generally offered in other cryptographic library implementations where cryptographic parameters are partially initialised.

- Offers unprecedented ease of use in a safe way. It is easy to put the group elements in a vector while ensuring that all vector elements are within the same type. It is easy to ensure that all inputs to a function are within the same type while easily allowing the swap of any cryptographic group. 

- It offers systematic safety by putting assertions within initialisation. This allows generic implementation of cryptographic libraries without needing to know about the underlying group; see, for instance, the DSA example, which is also ECDSA according to the FIPS standard.

- Memory safety via garbage collection of Julia's language ensures that bad things won't happen on this front;

- Disclaimer: Group operations are susceptible to side-channel timing attacks as constant cryptography is not currently used; hence, it is up to the developer to obfuscate runtime with delays. Alternatively, make contributions for constant curve implementations like Montgommery curve Curve25519.

- Performance: Projective coordinates for elliptic curves are currently not implemented, making them much slower than the state-of-the-art. Also, binary field implementations would highly benefit from using low-level assembly instructions for practical use.

- Offers cryptographic groups interface for all [FIPS 186.4](https://csrc.nist.gov/pubs/fips/186-4/final) curves, including:

  - Prime curves: P_192 (prime192v1, secp192r1), P_224, P_256 (prime256v1, secp256r1), P_384, P_521
  - Koblitz curves: K_163, K_233, K_283, K_409, K_571
  - Binary curves: B_163, B_233, B_288, B_409, B_571

  Koblitz and Binary curves can be instantiated in polynomial or Gaussian normal basis fields. Koblitz curves are mapped to generic binary curves, so the current implementation does not offer performance benefits. 

- Offers an easy way to instantiate any Weierstrass curve parameters; Montgomery and Edwards's curves are also in the scope of this package;

- Offers modular prime group constructor that offers safety at the constructor level via `g^q = 1` check;

- It is easy to experiment with cryptographic groups with macro constructors like `g = @PGroup{p=21, q=11}(2)` or `g = @ECGroup{P_192}()`, which also are used for aliasing type signatures, making stack traces much more pleasant—an ideal setup for anyone who is just starting their cryptographic journey.

## Elliptic Curve safety

- Elliptic curve safety is implemented in `ECPoint` wrapper type that adds essential safety checks;
- Elliptic curve point is rejected at `ECPoint` constructor if `p^h == 0` or if it does not satisfy the corresponding curve equation (see `oncurve` method);
- According to https://safecurves.cr.yp.to/complete.html, "The standard Weierstrass addition formulas fail if Q happens to match -P" is checked; 
