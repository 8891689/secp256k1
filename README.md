# secp256k1 Elliptic Curve Cryptography mini Library (C Implementation)

## Introduction

This project implements the core algorithms of secp256k1 elliptic curve cryptography in C.  The code does not rely on any external libraries and directly implements all necessary cryptographic primitives.

This library primarily implements cryptographic operations on the secp256k1 elliptic curve, with the following key features:

Big Integer Arithmetic: Provides basic operations for big integers (BigInt), including addition, subtraction, comparison, initialization, copying, zeroing, getting a specific bit, etc. These are the foundation for cryptographic calculations.

Montgomery Multiplication: Implements the Montgomery multiplication algorithm (montgomery_mult), an efficient method for modular multiplication used to accelerate elliptic curve operations. The library includes initialization of Montgomery parameters (montgomery_init). A more readable but less efficient mul_mod is also present.

Elliptic Curve Point Operations:

Supports point operations in both Jacobian (ECPointJac) and affine (ECPoint) coordinates.

Implements point addition (add_point_jac), point doubling (double_point_jac), and scalar multiplication (scalar_multiply_jac) in Jacobian coordinates. Scalar multiplication is optimized.

Implements point addition (point_add) and point doubling (double_point) in affine coordinates.

Provides conversion from Jacobian to affine coordinates (jacobian_to_affine).

Modular Arithmetic: Provides operations such as modular addition (add_mod), modular subtraction (sub_mod), modular exponentiation (modexp), and modular inverse (mod_inverse). There are also two versions of mod_generic and mul_mod.

Utility Functions:

print_bigint: Prints a big integer.

bigint_to_hex: Converts a big integer to a hexadecimal string.

hex_to_bigint: Converts a hexadecimal string to a big integer.

point_to_compressed_hex: Converts an elliptic curve point to a compressed hexadecimal string.

point_to_uncompressed_hex: Converts an elliptic curve point to an uncompressed hexadecimal string.

Summary:

This library is a relatively low-level cryptographic library for the secp256k1 curve. It does not provide high-level functionality like ECDSA signing/verification, but rather focuses on providing efficient, optimized elliptic curve point operations and modular arithmetic, serving as a foundation for building higher-level cryptographic applications (such as signing, key exchange, etc.). The library makes extensive use of manually optimized code (e.g., loop unrolling, pointer arithmetic) for better performance. It's important to notice that the Montgomery multiplication implementation is specifically optimized with respect to the secp256k1 curve parameters.

9-Word Extended Operations

- In modular multiplication, intermediate results (512 bits) are represented by an array of 9 `uint32_t`.
- Implements multiplication of a big integer by a constant, left shift by one word (32 bits), 9-word addition, and conversion from a 9-word array to `BigInt`. These functions are used in the optimization of modular multiplication.

Compilation

```bash
gcc -O3 -o test test.c secp256k1.c
```
```
./test
[Performance] 40000 iterations took: 1.012s (39542.0 ops/s)

Verification of the last 5 public keys:
Private key 39996:
  Hex: 0000000000000000000000000000000000000000000000000000000000009c3c
  Public Key: 03955b8ff94a2c9ae98acaa20300e1619619e7d86701888f1c440e9488212f261e
Private key 39997:
  Hex: 0000000000000000000000000000000000000000000000000000000000009c3d
  Public Key: 02ab6ab9880a79245f061a2458d8f1d61a3881ac9e38f9d3087d4853304957da74
Private key 39998:
  Hex: 0000000000000000000000000000000000000000000000000000000000009c3e
  Public Key: 03edbcc7280f833fb218547504958499de7fd593ebc7f6995ebbde624f4d6d0b47
Private key 39999:
  Hex: 0000000000000000000000000000000000000000000000000000000000009c3f
  Public Key: 02ccf949e4a9236746342d7856b199b6c5f0f9fa09202d4a2ceab9bfe6d9a71ff1
Private key 40000:
  Hex: 0000000000000000000000000000000000000000000000000000000000009c40
  Public Key: 0363cefbed3060c05aff68ce50aee6b8ea31f1ba496e23a3d43580c76795135814

```
