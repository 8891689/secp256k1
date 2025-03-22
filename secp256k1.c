/*Author: 8891689
  Assist in creation ：ChatGPT 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>

#include "secp256k1.h"

// ------------------------------
// BigInt 基础运算
// ------------------------------

void init_bigint(BigInt *x, uint32_t val) {
    memset(x->data, 0, sizeof(x->data));
    x->data[0] = val;
}

void bigint_init(BigInt *n) {
    memset(n->data, 0, sizeof(n->data));
}

void bigint_copy(const BigInt *src, BigInt *dest) {
    memcpy(dest->data, src->data, sizeof(src->data));
}

int compare_bigint(const BigInt *a, const BigInt *b) {
    for (int i = BIGINT_SIZE - 1; i >= 0; i--) {
        if (a->data[i] != b->data[i])
            return (a->data[i] > b->data[i]) ? 1 : -1;
    }
    return 0;
}

bool is_zero(const BigInt *a) {
    for (int i = 0; i < BIGINT_SIZE; i++) {
        if (a->data[i])
            return false;
    }
    return true;
}

int get_bit(const BigInt *a, int i) {
    return (a->data[i >> 5] >> (i & 31)) & 1;
}

// ------------------------------
// 低级加减运算
// ------------------------------

void ptx_u256Add(BigInt *res, const BigInt *a, const BigInt *b) {
    uint32_t carry = 0;
    for (int i = 0; i < BIGINT_SIZE; i++) {
        uint64_t sum = (uint64_t)a->data[i] + b->data[i] + carry;
        res->data[i] = (uint32_t)sum;
        carry = (uint32_t)(sum >> 32);
    }
}

void ptx_u256Sub(BigInt *res, const BigInt *a, const BigInt *b) {
    uint32_t borrow = 0;
    for (int i = 0; i < BIGINT_SIZE; i++) {
        uint64_t diff = (uint64_t)a->data[i] - b->data[i] - borrow;
        res->data[i] = (uint32_t)diff;
        borrow = (a->data[i] < b->data[i] + borrow) ? 1 : 0;
    }
}

// ------------------------------
// 9字扩展运算（用于模乘中间结果折叠）
// ------------------------------

void multiply_bigint_by_const(const BigInt *a, uint32_t c, uint32_t result[9]) {
    uint64_t carry = 0;
    for (int i = 0; i < BIGINT_SIZE; i++) {
        uint64_t prod = (uint64_t)a->data[i] * c + carry;
        result[i] = (uint32_t)prod;
        carry = prod >> 32;
    }
    result[BIGINT_SIZE] = (uint32_t)carry;
}

void shift_left_word(const BigInt *a, uint32_t result[9]) {
    result[0] = 0;
    memcpy(&result[1], a->data, BIGINT_SIZE * sizeof(uint32_t));
}

void add_9word(uint32_t r[9], const uint32_t addend[9]) {
    uint64_t carry = 0;
    for (int i = 0; i < 9; i++) {
        uint64_t sum = (uint64_t)r[i] + addend[i] + carry;
        r[i] = (uint32_t)sum;
        carry = sum >> 32;
    }
}

void convert_9word_to_bigint(const uint32_t r[9], BigInt *res) {
    memcpy(res->data, r, BIGINT_SIZE * sizeof(uint32_t));
}

// ------------------------------
// 模乘及模约简
// ------------------------------

void mul_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p) {
    uint32_t prod[16] = {0};
    for (int i = 0; i < BIGINT_SIZE; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < BIGINT_SIZE; j++) {
            uint64_t tmp = (uint64_t)prod[i+j] + (uint64_t)a->data[i] * b->data[j] + carry;
            prod[i+j] = (uint32_t)tmp;
            carry = tmp >> 32;
        }
        prod[i+BIGINT_SIZE] += (uint32_t)carry;
    }
    BigInt L, H;
    for (int i = 0; i < BIGINT_SIZE; i++) {
        L.data[i] = prod[i];
        H.data[i] = prod[i+BIGINT_SIZE];
    }
    uint32_t Rext[9] = {0};
    memcpy(Rext, L.data, BIGINT_SIZE * sizeof(uint32_t));
    Rext[BIGINT_SIZE] = 0;
    uint32_t H977[9] = {0};
    multiply_bigint_by_const(&H, 977, H977);
    add_9word(Rext, H977);
    uint32_t Hshift[9] = {0};
    shift_left_word(&H, Hshift);
    add_9word(Rext, Hshift);
    if (Rext[BIGINT_SIZE]) {
        uint32_t extra[9] = {0};
        BigInt extraBI;
        init_bigint(&extraBI, Rext[BIGINT_SIZE]);
        Rext[BIGINT_SIZE] = 0;
        uint32_t extra977[9] = {0}, extraShift[9] = {0};
        multiply_bigint_by_const(&extraBI, 977, extra977);
        shift_left_word(&extraBI, extraShift);
        memcpy(extra, extra977, 9 * sizeof(uint32_t));
        add_9word(extra, extraShift);
        add_9word(Rext, extra);
    }
    BigInt R_temp;
    convert_9word_to_bigint(Rext, &R_temp);
    if (Rext[BIGINT_SIZE] || compare_bigint(&R_temp, p) >= 0) {
        ptx_u256Sub(&R_temp, &R_temp, p);
        if (compare_bigint(&R_temp, p) >= 0)
            ptx_u256Sub(&R_temp, &R_temp, p);
    }
    bigint_copy(&R_temp, res);
}

void efficient_mod(BigInt *r, const BigInt *a, const BigInt *p) {
    bigint_copy(a, r);
    if (compare_bigint(r, p) >= 0) {
         BigInt temp;
         ptx_u256Sub(&temp, r, p);
         if (compare_bigint(&temp, p) >= 0)
              ptx_u256Sub(&temp, &temp, p);
         bigint_copy(&temp, r);
    }
}

void mod_generic(BigInt *r, const BigInt *a, const BigInt *p) {
    efficient_mod(r, a, p);
}

void sub_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p) {
    BigInt temp;
    if (compare_bigint(a, b) < 0) {
         BigInt sum;
         ptx_u256Add(&sum, a, p);
         ptx_u256Sub(&temp, &sum, b);
    } else {
         ptx_u256Sub(&temp, a, b);
    }
    mod_generic(&temp, &temp, p);
    bigint_copy(&temp, res);
}

void add_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p) {
    BigInt temp;
    ptx_u256Add(&temp, a, b);
    mod_generic(&temp, &temp, p);
    bigint_copy(&temp, res);
}

// ------------------------------
// 模幂与模逆
// ------------------------------

void modexp(BigInt *res, const BigInt *base, const BigInt *exp, const BigInt *p) {
    BigInt result;
    init_bigint(&result, 1);
    BigInt b;
    bigint_copy(base, &b);
    for (int i = 0; i < 256; i++) {
         if (get_bit(exp, i)) {
              BigInt temp;
              mul_mod(&temp, &result, &b, p);
              bigint_copy(&temp, &result);
         }
         BigInt temp;
         mul_mod(&temp, &b, &b, p);
         bigint_copy(&temp, &b);
    }
    bigint_copy(&result, res);
}

void mod_inverse(BigInt *res, const BigInt *a, const BigInt *p) {
    // 根据费马小定理：a^(p-2) mod p 为 a 的逆元
    BigInt p_minus_2;
    bigint_copy(p, &p_minus_2);
    BigInt two;
    init_bigint(&two, 2);
    BigInt temp;
    ptx_u256Sub(&temp, &p_minus_2, &two);
    bigint_copy(&temp, &p_minus_2);
    modexp(res, a, &p_minus_2, p);
}

// ------------------------------
// BigInt 与字节、整数转换
// ------------------------------

void bigint_from_bytes(BigInt *n, const uint8_t *bytes, size_t len) {
    bigint_init(n);
    size_t bytes_to_copy = (len > 32) ? 32 : len;
    for (size_t i = 0; i < bytes_to_copy; i++) {
        n->data[i / 4] |= ((uint32_t)bytes[bytes_to_copy - 1 - i]) << ((i % 4) * 8);
    }
}

void bigint_from_int(int num, BigInt *n) {
    bigint_init(n);
    if (num >= 0) {
        n->data[0] = (uint32_t)num;
    } else {
        n->data[0] = (uint32_t)(-num);
        for (int i = 0; i < BIGINT_SIZE; i++) {
            n->data[i] = ~n->data[i];
        }
        BigInt one;
        bigint_init(&one);
        one.data[0] = 1;
        BigInt temp;
        bigint_add(n, &one, &temp, NULL);
        bigint_copy(&temp, n);
    }
}

void bigint_divide_by_const(const BigInt *a, uint32_t c, BigInt *quotient) {
    uint64_t r = 0;
    for (int i = BIGINT_SIZE - 1; i >= 0; i--) {
        r = (r << 32) | a->data[i];
        quotient->data[i] = (uint32_t)(r / c);
        r = r % c;
    }
}

// bigint_add: 实现 a + b，若 mod 非空则做一次模约简
void bigint_add(const BigInt *a, const BigInt *b, BigInt *res, const BigInt *mod) {
    BigInt sum;
    ptx_u256Add(&sum, a, b);
    if (mod) {
        if (compare_bigint(&sum, mod) >= 0) {
            BigInt temp;
            ptx_u256Sub(&temp, &sum, mod);
            bigint_copy(&temp, res);
        } else {
            bigint_copy(&sum, res);
        }
    } else {
        bigint_copy(&sum, res);
    }
}

// bigint_mul: 若 mod 非空则使用 mul_mod，否则仅计算低256位结果
void bigint_mul(const BigInt *a, const BigInt *b, BigInt *res, const BigInt *mod) {
    if (mod) {
        mul_mod(res, a, b, mod);
    } else {
        uint32_t prod[16] = {0};
        for (int i = 0; i < BIGINT_SIZE; i++) {
            uint64_t carry = 0;
            for (int j = 0; j < BIGINT_SIZE; j++) {
                uint64_t tmp = (uint64_t)prod[i+j] + (uint64_t)a->data[i] * b->data[j] + carry;
                prod[i+j] = (uint32_t)tmp;
                carry = tmp >> 32;
            }
            prod[i+BIGINT_SIZE] += (uint32_t)carry;
        }
        for (int i = 0; i < BIGINT_SIZE; i++) {
            res->data[i] = prod[i];
        }
    }
}

// ------------------------------
// 模平方根 sqrt_mod_p
// ------------------------------
// 适用于 secp256k1 素数 p (p % 4 == 3)，使用公式：sqrt(a) = a^((p+1)/4) mod p
void sqrt_mod_p(BigInt *res, const BigInt *a, const BigInt *p) {
    BigInt p_plus_1, exp;
    bigint_copy(p, &p_plus_1);  // p_plus_1 = p
    BigInt one;
    bigint_init(&one);
    one.data[0] = 1;
    BigInt tmp;
    ptx_u256Add(&tmp, &p_plus_1, &one);  // tmp = p + 1
    bigint_copy(&tmp, &p_plus_1);
    // exp = (p_plus_1) / 4
    bigint_divide_by_const(&p_plus_1, 4, &exp);
    // res = a^exp mod p
    modexp(res, a, &exp, p);
}

// ------------------------------
// ECC 点运算
// ------------------------------

void point_set_infinity(ECPoint *P) {
    P->infinity = true;
}

void point_copy(ECPoint *dest, const ECPoint *src) {
    bigint_copy(&src->x, &dest->x);
    bigint_copy(&src->y, &dest->y);
    dest->infinity = src->infinity;
}

void point_add(ECPoint *R, const ECPoint *P, const ECPoint *Q, const BigInt *p) {
    if (P->infinity) { point_copy(R, Q); return; }
    if (Q->infinity) { point_copy(R, P); return; }
    BigInt diffY, diffX, inv_diffX, lambda, lambda2, temp;
    sub_mod(&diffY, &Q->y, &P->y, p);
    sub_mod(&diffX, &Q->x, &P->x, p);
    mod_inverse(&inv_diffX, &diffX, p);
    mul_mod(&lambda, &diffY, &inv_diffX, p);
    mul_mod(&lambda2, &lambda, &lambda, p);
    // x_R = lambda^2 - x_P - x_Q
    sub_mod(&temp, &lambda2, &P->x, p);
    sub_mod(&R->x, &temp, &Q->x, p);
    // y_R = lambda * (x_P - x_R) - y_P
    sub_mod(&temp, &P->x, &R->x, p);
    mul_mod(&R->y, &lambda, &temp, p);
    sub_mod(&R->y, &R->y, &P->y, p);
    R->infinity = false;
}

void double_point(ECPoint *R, const ECPoint *P, const BigInt *p) {
    if (P->infinity || is_zero(&P->y)) {
         point_set_infinity(R);
         return;
    }
    BigInt x2, numerator, denominator, inv_den, lambda, lambda2, two, two_x;
    mul_mod(&x2, &P->x, &P->x, p);
    BigInt three; 
    init_bigint(&three, 3);
    mul_mod(&numerator, &three, &x2, p);
    init_bigint(&two, 2);
    mul_mod(&denominator, &two, &P->y, p);
    mod_inverse(&inv_den, &denominator, p);
    mul_mod(&lambda, &numerator, &inv_den, p);
    mul_mod(&lambda2, &lambda, &lambda, p);
    mul_mod(&two_x, &two, &P->x, p);
    sub_mod(&R->x, &lambda2, &two_x, p);
    sub_mod(&numerator, &P->x, &R->x, p);
    mul_mod(&R->y, &lambda, &numerator, p);
    sub_mod(&R->y, &R->y, &P->y, p);
    R->infinity = false;
}

void scalar_multiply(const BigInt *d, const ECPoint *G, const BigInt *p, ECPoint *result) {
    ECPoint R;
    point_set_infinity(&R);  // R = O
    for (int i = 255; i >= 0; i--) {
        ECPoint temp;
        double_point(&temp, &R, p);
        point_copy(&R, &temp);
        if (get_bit(d, i)) {
            ECPoint temp2;
            point_add(&temp2, &R, G, p);
            point_copy(&R, &temp2);
        }
    }
    point_copy(result, &R);
}

// ------------------------------
// 辅助工具函数
// ------------------------------

void print_bigint(const BigInt *b) {
    for (int i = BIGINT_SIZE - 1; i >= 0; i--) {
        printf("%08x", b->data[i]);
    }
    printf("\n");
}

void hex_to_bigint(const char *hex, BigInt *b) {
    memset(b->data, 0, sizeof(b->data));
    int len = (int)strlen(hex);
    int j = 0;
    for (int i = len; i > 0; i -= 8) {
        int start = i - 8;
        if (start < 0) start = 0;
        char temp[9] = {0};
        strncpy(temp, hex + start, 8);
        b->data[j++] = (uint32_t)strtoul(temp, NULL, 16);
    }
}

void bigint_to_hex(const BigInt *num, char *hex_string) {
    int len = 0;
    for (int i = BIGINT_SIZE - 1; i >= 0; i--) {
        len += sprintf(hex_string + len, "%08x", num->data[i]);
    }
}

void point_to_compressed_hex(const ECPoint *P, char *hex_string) {
    if (P->infinity) {
        strcpy(hex_string, "00");
        return;
    }
    char x_hex[65];
    bigint_to_hex(&P->x, x_hex);
    for (int i = 0; x_hex[i]; i++) {
        x_hex[i] = tolower(x_hex[i]);
    }
    if ((P->y.data[0] & 1) == 0)
        sprintf(hex_string, "02%s", x_hex);
    else
        sprintf(hex_string, "03%s", x_hex);
}

void point_to_uncompressed_hex(const ECPoint *P, char *hex_string) {
    if (P->infinity) {
        strcpy(hex_string, "00");
        return;
    }
    char x_hex[65], y_hex[65];
    bigint_to_hex(&P->x, x_hex);
    bigint_to_hex(&P->y, y_hex);
    for (int i = 0; x_hex[i]; i++) {
        x_hex[i] = tolower(x_hex[i]);
    }
    for (int i = 0; y_hex[i]; i++) {
        y_hex[i] = tolower(y_hex[i]);
    }
    sprintf(hex_string, "04%s%s", x_hex, y_hex);
}
