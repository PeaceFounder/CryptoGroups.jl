# CryptoGroups.jl

[![codecov](https://codecov.io/gh/PeaceFounder/CryptoGroups.jl/graph/badge.svg?token=G9HT9VSV4T)](https://codecov.io/gh/PeaceFounder/CryptoGroups.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PeaceFounder.github.io/CryptoGroups.jl/dev)

CryptoGroups is a Julia package that provides a versatile and type-safe implementation of cryptographic groups. It offers a unified interface for working with various types of groups, including modular prime groups and elliptic curves over prime and binary fields. Suitable for both academic and aspiring production environments, CryptoGroups has been rigorously tested and refined through its implementation in real-world projects such as [CryptoSignatures](https://github.com/PeaceFounder/CryptoSignatures.jl) and [ShuffleProofs](https://github.com/PeaceFounder/ShuffleProofs.jl), ensuring a robust and versatile cryptographic toolkit.

## Key Features

1. **Group-Agnostic API**: CryptoGroups allows for polymorphic implementations of cryptographic standards, such as the [FIPS 186-4](https://csrc.nist.gov/pubs/fips/186-4/final) digital signature algorithm, that work seamlessly across different group types.

2. **Type Safety**: The package leverages Julia's powerful type system, encoding relevant group parameters as type parameter values. This design minimizes the need for runtime assertions and ensures efficient memory use.

3. **Comprehensive Group Support**: CryptoGroups supports a wide range of cryptographic groups, including:
   - Modular prime groups
   - Elliptic curves over prime fields
   - Elliptic curves over binary fields
   - All curves specified in [FIPS 186-4](https://csrc.nist.gov/pubs/fips/186-4/final), including Weierstrass curves, Koblitz curves, and Pseudorandom curves:

   | Curve Type | Curves |
   |------------|--------|
   | Weierstrass Curves | `P_192` (`prime192v1`, `secp192r1`), `P_224`, `P_256` (`prime256v1`, `secp256r1`), `P_384`, `P_521` |
   | Koblitz Curves | `K_163`, `K_233`, `K_283`, `K_409`, `K_571` |
   | Pseudorandom Curves | `B_163`, `B_233`, `B_288`, `B_409`, `B_571` |

   Note: Koblitz and Pseudorandom curves can be instantiated in polynomial or Gaussian normal basis (default) fields. Currently, Koblitz curves are mapped to generic binary curves, so there are no performance benefits specific to Koblitz curves in the current implementation. In future support for Montgomery and Edwards curves as specified in [FIPS 186-5](https://csrc.nist.gov/pubs/fips/186-5/final) may be added. 

4. **Safety-First Approach**: Essential security checks, such as verifying point membership on elliptic curves or cofactor validation, are integrated into the constructors. Special cases in arithmetic are treated properly.

5. **Easy-to-Use Interface**: CryptoGroups provides macro constructors for quick experimentation and type aliasing, making it ideal for both beginners and experienced cryptographers.

## Main Components

### Group Types

1. `Group`: An abstract type representing a cyclic group.
2. `PGroup`: Represents modular prime groups (instantiated with `@PGRoup`).
3. `ECGroup`: Represents elliptic curve groups (instantiated with `@ECGRoup`).

### Key Methods

- `order`: Get the order of the group
- `octet`: Convert group elements to octet representation
- `value`: Convert group elements to their simplest representation
- `*`, `/`, `^`: Perform group operations (multiplication, division, exponentiation)
- `inv`: Compute the inverse of a group element
- `one`: Construct the identity element of the group

### Utility Functions

- `concretize_type`: Construct concrete subtypes of abstract group types
- `spec`: Get or construct group specifications
- `iscompressable`: Check if a group element can be compressed in octet representation

## Usage Examples

1. Creating a modular prime group:
   ```julia
   G = @PGroup{p = 23, q = 11}
   g = G(2)
   ```

2. Creating an elliptic curve group:
   ```julia
   G = @ECGroup{P_192}
   g = G() # uses generator from specification
   ```

3. Performing group operations:
   ```julia
   g^3 * g^5 / g^2 == (g^3)^2 == g^6
   g^(order(G) - 1) * g == one(G)
   ```
   
4. Serializing and deserializing:
   ```julia
   g == G(octet(g)) == G(value(g))
   ```

## Safety Guarantees

While no cryptographic system can guarantee absolute security, CryptoGroups implements the following safety features:

- Group element arithmetics is possible only with the same types of groups and throws `MethodError` when that is violated. For instance, `@ECGroup{P_192}() * @ECGroup{P_256}()` throws an error;
- Group elements are validated during construction, throwing `ArgumentError` for invalid inputs;
- Modular prime group elements are checked to belong in prime group via $g^q = 1$ or with an efficient `jacobi(g, p) = 1` for quadratic residue groups;
- Elliptic curve points are checked for curve equation satisfaction and cofactor validation;
- The package implements checks to prevent issues with special cases in point addition formulas;
- Exponentiation with $k~ {\rm mod} ~q = 0$ shows warning or throws an error in a strict mode.

## Limitations and Future Work

CryptoGroups prioritizes type safety and flexibility in its implementation of cryptographic groups. However, this focus comes with some trade-offs in terms of performance, which in turn can have implications for security in certain scenarios.

The current implementation of CryptoGroups has several areas where performance optimizations are yet to be implemented:

- The package doesn't use projective coordinates for elliptic curve arithmetics;
- Lacks special treatment for Koblitz curves;
- Doesn't implement Mersenne primes when available over generic prime fields;
- Binary field operations, the current implementation is suboptimal and doesn't take advantage of hardware-provided carryless operations.

These limitations result in significantly slower performance compared to state-of-the-art implementations. Preliminary estimates suggest that operations on prime curves in CryptoGroups are about 100 times slower than optimized libraries like OpenSSL, while binary curves may be up to 1000 times slower.

These performance limitations can potentially lead to security risks in certain contexts. The most immediate concern is the increased vulnerability to Denial of Service (DoS) attacks. Since computations take longer to complete, systems using CryptoGroups for cryptographic operations could be more easily overwhelmed by a flood of requests. Additionally, the current arithmetic operations are not side-channel resistant, which could pose a risk if an adversary is able to monitor the machine performing the group operations.

For scenarios where only remote adversaries are a concern, the side-channel risk can be mitigated via time-padding. Set a fixed computation time based on benchmarks, add sleep periods to reach this time, and reject results if the computation exceeds the set time due to system overload. However, this approach doesn't address the fundamental performance limitations and may not be suitable for all use cases.

Looking ahead, there are several areas where CryptoGroups could be improved to address these performance and security concerns. Future work should focus on implementing projective coordinates for elliptic curves, optimizing binary field operations using hardware-specific instructions, and developing constant-time arithmetic operations to improve side-channel resistance. Additionally, conducting thorough security audits and expanding test coverage will be crucial to identify and address any overlooked vulnerabilities.

Despite these limitations, CryptoGroups provides a solid foundation for cryptographic operations in Julia. Its type-safe design and comprehensive group support make it an excellent tool for educational purposes, allowing users to experiment with and learn about various cryptographic groups. For production use, while CryptoGroups can be deployed, users should carefully consider the performance implications and implement necessary security measures. With appropriate caution and mitigation strategies, CryptoGroups can be a valuable component in cryptographic applications, particularly in scenarios where its strengths in type safety and flexibility outweigh the performance limitations.

# References

- [elliptic-curve](https://github.com/sdiehl/elliptic-curve#readme) library in Haskell that shares similar goals regarding type safety
- [RFC2409](https://tools.ietf.org/html/rfc2409#section-6.2) and [RFC5114](https://tools.ietf.org/html/rfc5114#section-2.1) for modular prime group constants
- [SafeCurves](https://safecurves.cr.yp.to/complete.html) on addition checks for Weierstrass curves
- [FIPS 186-4](https://csrc.nist.gov/pubs/fips/186-4/final) and [FIPS 186-5](https://csrc.nist.gov/pubs/fips/186-5/final)
- [NIST SP 800-186](https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-186-draft.pdf)
- [ANSI X9.142](https://webstore.ansi.org/preview-pages/ASCX9/preview_ANSI+X9.142-2020.pdf) and in unpaywalled form [here](https://www.cs.miami.edu/home/burt/learning/Csc609.142/ecdsa-cert.pdf)
- [CryptoSignatures.jl](https://github.com/PeaceFounder/CryptoSignatures.jl) FIPS 186-4 digital signature algorithm implemetation
- [ShuffleProofs.jl](https://github.com/PeaceFounder/ShuffleProofs.jl) Verificatum compatable ElGamal proof of shuffle implementation
