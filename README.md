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
Generated 20000 public keys in 1.3837 seconds
Last 5 public keys and corresponding private keys:
Private Key 19996: 0000000000000000000000000000000000000000000000000000000000004e1c
Compressed Public Key: 0298244597f39daac5ca35a0769aeff54c5caa37586e260adb0596c65962f238c6
Uncompressed Public Key: 0498244597f39daac5ca35a0769aeff54c5caa37586e260adb0596c65962f238c6c698a22e103781276f04dcb519afd310f42e888b00031ee0a5565a6ecb988732

Private Key 19997: 0000000000000000000000000000000000000000000000000000000000004e1d
Compressed Public Key: 0299cb1469e9aa44880618c30c6eb8dad7717494008695225b592a3fe4f59cbc24
Uncompressed Public Key: 0499cb1469e9aa44880618c30c6eb8dad7717494008695225b592a3fe4f59cbc2405e404dc7091100fc29231e95ce2f3c5e228b9abda3af6e269ff71ac9a9040c6

Private Key 19998: 0000000000000000000000000000000000000000000000000000000000004e1e
Compressed Public Key: 0349421e70e35a26b188156bd1408ad7524db900f8085fbc5b5fbc91c1602872fd
Uncompressed Public Key: 0449421e70e35a26b188156bd1408ad7524db900f8085fbc5b5fbc91c1602872fdd32acfff1c1627b5eaed27246834e352913dc9d6d280f8890336edf2c5b62517

Private Key 19999: 0000000000000000000000000000000000000000000000000000000000004e1f
Compressed Public Key: 025af87006e0e4f30bb44506306e789a2b7ce6862fa9a419e4badd5450cccb074d
Uncompressed Public Key: 045af87006e0e4f30bb44506306e789a2b7ce6862fa9a419e4badd5450cccb074df2dee17cde9684e1d412ffd6ce54662fcfac45221504ec97904d52eb19ff69de

Private Key 20000: 0000000000000000000000000000000000000000000000000000000000004e20
Compressed Public Key: 021a03aba6c1d66b0adff5f523b05ae59226b75a3c89c5755728d4278b4d02dec0
Uncompressed Public Key: 041a03aba6c1d66b0adff5f523b05ae59226b75a3c89c5755728d4278b4d02dec0a7a83297eda03a1795ae0917fbd47ef29875e655eb0083e8421f27d1f5c06f9a

```
