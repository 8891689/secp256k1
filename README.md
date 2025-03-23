# secp256k1 Elliptic Curve Cryptography mini Library (C Implementation)

## Introduction

This project implements the core algorithms of secp256k1 elliptic curve cryptography in C.  The code does not rely on any external libraries and directly implements all necessary cryptographic primitives.

## Implementation Details

### Big Integer Arithmetic

- 256-bit unsigned integers are represented using an array of 8 `uint32_t` (`BigInt` struct).
- Implements addition, subtraction, multiplication, comparison, and bitwise operations (getting a specific bit) for big integers.
- Addition and subtraction simulate manual calculation, handling carries and borrows.
- Big integer multiplication uses the basic schoolbook multiplication algorithm, producing a 512-bit result.
- Comparison compares words from most significant to least significant.
- Bitwise operations are implemented using shifts and masks.

### Modular Arithmetic

- **Modular multiplication (`mul_mod`) uses an optimized implementation**:
    1.  First, the product of two 256-bit big integers (512-bit) is calculated.
    2.  The 512-bit result is split into the high 256 bits (H) and the low 256 bits (L).
    3.  Fast reduction is performed using the special form of the secp256k1 curve modulus `p` (p = 2^256 - 2^32 - 977):
        - Multiply H by 977, the result is H977.
        - Shift H left by one uint32_t word (32 bits), the result is Hshift.
        - Add L, H977, and Hshift (9 `uint32_t` words long), the result is Rext.
        - If Rext has a ninth bit (overflow bit), perform extra processing:
            - Multiply the overflow by 977 to get extra977, and shift left to get extraShift.
            - Add Rext, extra977, and extraShift.
        - Convert Rext to `BigInt`.
    4.  If the result is still greater than or equal to `p`, perform at most two subtractions.
- Modular addition (`add_mod`) and modular subtraction (`sub_mod`) perform one or two modular subtractions based on big integer addition and subtraction.
- Modular inverse (`mod_inverse`) is implemented using Fermat's Little Theorem and modular exponentiation (a<sup>p-2</sup> mod p).
- Modular exponentiation (`modexp`) uses the "left-to-right" square-and-multiply algorithm.
- Optimized modular reduction function (`efficient_mod`) and generic modular reduction function (`mod_generic`) are provided. In the optimized version, if the comparison result is greater than or equal to p, subtract p at most twice.

### Elliptic Curve Point Operations

- The `ECPoint` struct represents a point on the elliptic curve, containing x and y coordinates and a flag indicating whether it is the point at infinity.
- Point addition (`point_add`) and point doubling (`double_point`) operations are implemented according to standard formulas, handling the point at infinity.
- Scalar multiplication (`scalar_multiply`) uses the classic "double-and-add" algorithm:
    - Initialize the result point R to infinity.
    - Iterate through the bits of the scalar from most significant to least significant:
        - Double R (`double_point`).
        - If the current bit is 1, add P to R (`point_add`).
    - The final R is the result.

### Utility Functions

- `print_bigint`: Prints a `BigInt` in hexadecimal format.
- `hex_to_bigint`: Converts a hexadecimal string to a `BigInt`.
- `bigint_to_hex`: Converts a `BigInt` to a hexadecimal string.

### 9-Word Extended Operations

- In modular multiplication, intermediate results (512 bits) are represented by an array of 9 `uint32_t`.
- Implements multiplication of a big integer by a constant, left shift by one word (32 bits), 9-word addition, and conversion from a 9-word array to `BigInt`. These functions are used in the optimization of modular multiplication.

## Compilation

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
