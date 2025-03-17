/*Author: 8891689/ChatGPT */
#include "secp256k1.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ------------------------------
// 大整数基础运算实现
// ------------------------------

void init_bigint(BigInt *x, uint32_t val) {
    memset(x->data, 0, sizeof(x->data));
    x->data[0] = val;
}

void copy_bigint(BigInt *dest, const BigInt *src) {
    memcpy(dest->data, src->data, sizeof(src->data));
}

int compare_bigint(const BigInt *a, const BigInt *b) {
    for (int i = 7; i >= 0; i--) {
        if (a->data[i] != b->data[i])
            return (a->data[i] > b->data[i]) ? 1 : -1;
    }
    return 0;
}

bool is_zero(const BigInt *a) {
    for (int i = 0; i < 8; i++) {
        if (a->data[i])
            return false;
    }
    return true;
}

int get_bit(const BigInt *a, int i) {
    // i>>5 代替除以32，i&31 代替取余
    return (a->data[i >> 5] >> (i & 31)) & 1;
}

void ptx_u256Add(BigInt *res, const BigInt *a, const BigInt *b) {
    uint32_t carry = 0;
    for (int i = 0; i < 8; i++) {
        uint64_t sum = (uint64_t)a->data[i] + b->data[i] + carry;
        res->data[i] = (uint32_t)sum;
        carry = (uint32_t)(sum >> 32);
    }
}

void ptx_u256Sub(BigInt *res, const BigInt *a, const BigInt *b) {
    uint32_t borrow = 0;
    for (int i = 0; i < 8; i++) {
        uint64_t diff = (uint64_t)a->data[i] - b->data[i] - borrow;
        res->data[i] = (uint32_t)diff;
        borrow = (a->data[i] < b->data[i] + borrow) ? 1 : 0;
    }
}

// ------------------------------
// 9字扩展运算实现（用于模乘中中间结果的折叠）
// ------------------------------

void multiply_bigint_by_const(const BigInt *a, uint32_t c, uint32_t result[9]) {
    uint64_t carry = 0;
    for (int i = 0; i < 8; i++) {
        uint64_t prod = (uint64_t)a->data[i] * c + carry;
        result[i] = (uint32_t)prod;
        carry = prod >> 32;
    }
    result[8] = (uint32_t)carry;
}

void shift_left_word(const BigInt *a, uint32_t result[9]) {
    result[0] = 0;
    memcpy(&result[1], a->data, 8 * sizeof(uint32_t));
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
    memcpy(res->data, r, 8 * sizeof(uint32_t));
}

// ------------------------------
// 模乘及模约简相关实现
// ------------------------------

void mul_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p) {
    uint32_t prod[16] = {0};
    // 计算 256 位大整数相乘得到 512 位结果
    for (int i = 0; i < 8; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 8; j++) {
            uint64_t tmp = (uint64_t)prod[i+j] + (uint64_t)a->data[i] * b->data[j] + carry;
            prod[i+j] = (uint32_t)tmp;
            carry = tmp >> 32;
        }
        prod[i+8] += (uint32_t)carry;
    }
    BigInt L, H;
    for (int i = 0; i < 8; i++) {
        L.data[i] = prod[i];
        H.data[i] = prod[i+8];
    }
    uint32_t Rext[9] = {0};
    memcpy(Rext, L.data, 8 * sizeof(uint32_t));
    Rext[8] = 0;
    uint32_t H977[9] = {0};
    multiply_bigint_by_const(&H, 977, H977);
    add_9word(Rext, H977);
    uint32_t Hshift[9] = {0};
    shift_left_word(&H, Hshift);
    add_9word(Rext, Hshift);
    if (Rext[8]) {
        uint32_t extra[9] = {0};
        BigInt extraBI;
        init_bigint(&extraBI, Rext[8]);
        Rext[8] = 0;
        uint32_t extra977[9] = {0}, extraShift[9] = {0};
        multiply_bigint_by_const(&extraBI, 977, extra977);
        shift_left_word(&extraBI, extraShift);
        memcpy(extra, extra977, 9 * sizeof(uint32_t));
        add_9word(extra, extraShift);
        add_9word(Rext, extra);
    }
    BigInt R_temp;
    convert_9word_to_bigint(Rext, &R_temp);
    if (Rext[8] || compare_bigint(&R_temp, p) >= 0) {
        ptx_u256Sub(&R_temp, &R_temp, p);
        if (compare_bigint(&R_temp, p) >= 0)
            ptx_u256Sub(&R_temp, &R_temp, p);
    }
    copy_bigint(res, &R_temp);
}

void efficient_mod(BigInt *r, const BigInt *a, const BigInt *p) {
    copy_bigint(r, a);
    if (compare_bigint(r, p) >= 0) {
         BigInt temp;
         ptx_u256Sub(&temp, r, p);
         if (compare_bigint(&temp, p) >= 0)
              ptx_u256Sub(&temp, &temp, p);
         copy_bigint(r, &temp);
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
    mod_generic(res, &temp, p);
}

void add_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p) {
    BigInt temp;
    ptx_u256Add(&temp, a, b);
    mod_generic(res, &temp, p);
}

void modexp(BigInt *res, const BigInt *base, const BigInt *exp, const BigInt *p) {
    BigInt result;
    init_bigint(&result, 1);
    BigInt b;
    copy_bigint(&b, base);
    for (int i = 0; i < 256; i++) {
         if (get_bit(exp, i)) {
              BigInt temp;
              mul_mod(&temp, &result, &b, p);
              copy_bigint(&result, &temp);
         }
         BigInt temp;
         mul_mod(&temp, &b, &b, p);
         copy_bigint(&b, &temp);
    }
    copy_bigint(res, &result);
}

void mod_inverse(BigInt *res, const BigInt *a, const BigInt *p) {
    // 根据费马小定理：a^(p-2) mod p 为 a 的逆元
    BigInt p_minus_2;
    copy_bigint(&p_minus_2, p);
    BigInt two;
    init_bigint(&two, 2);
    BigInt temp;
    ptx_u256Sub(&temp, &p_minus_2, &two);
    copy_bigint(&p_minus_2, &temp);
    modexp(res, a, &p_minus_2, p);
}

// ------------------------------
// ECC 点运算实现
// ------------------------------

void point_set_infinity(ECPoint *P) {
    P->infinity = true;
}

void point_copy(ECPoint *dest, const ECPoint *src) {
    copy_bigint(&dest->x, &src->x);
    copy_bigint(&dest->y, &src->y);
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

// 使用传统双倍加法算法实现标量乘法
void scalar_multiply(const BigInt *d, const ECPoint *P, const BigInt *p, ECPoint *result) {
    ECPoint R;
    point_set_infinity(&R);  // R = O
    for (int i = 255; i >= 0; i--) {
        // R = 2R
        ECPoint temp;
        double_point(&temp, &R, p);
        point_copy(&R, &temp);
        if (get_bit(d, i)) {
            // R = R + P
            ECPoint temp2;
            point_add(&temp2, &R, P, p);
            point_copy(&R, &temp2);
        }
    }
    point_copy(result, &R);
}

// ------------------------------
// 辅助工具函数实现
// ------------------------------

void print_bigint(const BigInt *b) {
    for (int i = 7; i >= 0; i--) {
        printf("%08x", b->data[i]);
    }
    printf("\n");
}

void hex_to_bigint(const char* hex, BigInt *b) {
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
    for (int i = BIGINT_SIZE - 1; i >= 0; i--) { // 现在 BIGINT_SIZE 已定义
        len += sprintf(hex_string + len, "%08X", num->data[i]);
    }

    // 移除前导零
    int start = 0;
    while (hex_string[start] == '0' && hex_string[start + 1] != '\0') {
        start++;
    }
    if (start > 0) {
        memmove(hex_string, hex_string + start, strlen(hex_string) - start + 1);
    }
    if (hex_string[0] == '\0') { // 如果全为零
        strcpy(hex_string, "0");
    }
}

