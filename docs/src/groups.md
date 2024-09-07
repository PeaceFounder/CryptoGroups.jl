# Groups

The CryptoGroups package provides a main API for working with cryptographic groups in a group-agnostic manner, supporting modular prime groups, elliptic curves over prime fields, and binary fields. This flexibility enables polymorphic implementations of standards like FIPS 186-4 digital signature algorithm to work seamlessly across different group types. Essential security checks, such as verifying point membership on elliptic curves or cofactor validation, are integrated into the constructors.

CryptoGroups encodes relevant group parameters as type parameter values, leveraging Julia's powerful type system to minimize the need for runtime assertions and ensure efficient memory use. This design allows for vector-like collections to be used without subtyping AbstractArray, resulting in a leaner and more predictable codebase. This also necessitates compiling every group, and are up to the user to decide which groups are needed for their applications.

```@autodocs
Modules = [CryptoGroups]
Order = [:type, :function]
```

