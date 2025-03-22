# secp256k1 Elliptic Curve Cryptography Library (C Implementation)

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
./test
生成50000个公钥耗时: 3.5878秒
最后5个公钥坐标:
公钥49996 x: 0DB4754A9225CB3027B1327A80057B40286ACBF034FB4585B8B7734A1F204B57
公钥49996 y: 991BA0D69B3E3A9CF773F5FA59472D0B462C6D8964DB4CC34EF7DCE9565A5510
公钥49997 x: 42F21B577C238B1C22B2329193212D5971643FBBF95D8AA9381209B5505D850E
公钥49997 y: 2267BBD290783C9CF8743C883FF5491FFF2B20C2E49CC13DE543E583A18A4E1B
公钥49998 x: F4882926C56EEF04DD47DF0F565D72B84C76230415AE532D4B5F6FB9492428E8
公钥49998 y: 614D9B00487230B6A95F99CA79591C0F7489357CDA8D2B3139542E05FE78D981
公钥49999 x: 725AB6FC121A937084BE5AC4A4EBA0815026826B9AA0E5EF28753E2667181E94
公钥49999 y: 5524D38FCC64DD0A7E20649496BF18BCA2671187BB886ABF63C19DC4DFB844DD
公钥50000 x: AEB457CF4DCB83A294F187A06397BDCEF61F3AF5ABC5A866904F13B56436CD03
公钥50000 y: 3885DC959132FD17D761DA1C6B52EE841B924E12F05E01EACD18505BF96D874D
