# Groups

The CryptoGroups exposes main API to work with cryptographic groups in underlying group agnostic way be  it a modulus prime group, ellliptic curve over prime or binary fields. This for instance enables to implement FIPS 186.4 standard digital signature algorithm comningin modulus prime group and elliptic curve implementations. The underlying assertions like checking that point belongs to elliptic curve or it's cofactor but are put within the the constructor. On top of that cryptographic groups encodes all relevant parameters in type parameters allowing use type system instead of burdening runtime assertions and ensures efficient store within vector like collections without needing to subtype `AbstractArray` making codebase leaner and more prodeicatble. 

```@autodocs
Modules = [CryptoGroups]
Order = [:type, :function]
```

