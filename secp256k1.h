/*Author: 8891689
 * Assist in creation ：ChatGPT 
 */
#ifndef SECP256K1_H
#define SECP256K1_H

#include <stdint.h>
#include <stdbool.h>

// 定义大整数，每个大整数由 8 个 32 位无符号整数构成
#define BIGINT_WORDS 8

typedef struct {
    uint32_t data[BIGINT_WORDS];
} BigInt;

// 仿射坐标下的 ECC 点
typedef struct {
    BigInt x, y;
    bool infinity;
} ECPoint;

// 雅可比坐标下的 ECC 点
typedef struct {
    BigInt X, Y, Z;
    bool infinity;
} ECPointJac;

// ------------------------------
// 大整数运算函数原型
// ------------------------------
void init_bigint(BigInt *x, uint32_t val);
void copy_bigint(BigInt *dest, const BigInt *src);
int compare_bigint(const BigInt *a, const BigInt *b);
bool is_zero(const BigInt *a);
int get_bit(const BigInt *a, int i);
void ptx_u256Add(BigInt *res, const BigInt *a, const BigInt *b);
void ptx_u256Sub(BigInt *res, const BigInt *a, const BigInt *b);

// 扩展运算
void multiply_bigint_by_const(const BigInt *a, uint32_t c, uint32_t result[9]);
void shift_left_word(const BigInt *a, uint32_t result[9]);
void add_9word(uint32_t r[9], const uint32_t addend[9]);
void convert_9word_to_bigint(const uint32_t r[9], BigInt *res);

// 模乘及模约简
void mul_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p);
void efficient_mod(BigInt *r, const BigInt *a, const BigInt *p);
void mod_generic(BigInt *r, const BigInt *a, const BigInt *p);
void sub_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p);
void add_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p);
void modexp(BigInt *res, const BigInt *base, const BigInt *exp, const BigInt *p);
void mod_inverse(BigInt *res, const BigInt *a, const BigInt *p);

// ------------------------------
// 仿射坐标下 ECC 点运算函数原型
// ------------------------------
void point_set_infinity(ECPoint *P);
void point_copy(ECPoint *dest, const ECPoint *src);
void point_add(ECPoint *R, const ECPoint *P, const ECPoint *Q, const BigInt *p);
void double_point(ECPoint *R, const ECPoint *P, const BigInt *p);

// ------------------------------
// 雅可比坐标下 ECC 点运算函数原型
// ------------------------------
void point_set_infinity_jac(ECPointJac *P);
void point_copy_jac(ECPointJac *dest, const ECPointJac *src);
void double_point_jac(ECPointJac *R, const ECPointJac *P, const BigInt *p);
void add_point_jac(ECPointJac *R, const ECPointJac *P, const ECPointJac *Q, const BigInt *p);
void jacobian_to_affine(ECPoint *R, const ECPointJac *P, const BigInt *p);

// 内存管理
void free_bigint(BigInt *a);
void free_point_jac(ECPointJac *p);

// 雅可比坐标初始化
void init_point_jac(ECPointJac *p, bool infinity);

// 雅可比坐标复制
void copy_point_jac(ECPointJac *dest, const ECPointJac *src);

// 标量乘法核心函数
void scalar_multiply_jac(ECPointJac *result, 
                        const ECPointJac *point,
                        const BigInt *scalar,
                        const BigInt *p);

// --- 辅助工具 ---
void print_bigint(const BigInt *b);
void hex_to_bigint(const char *hex, BigInt *b);
void bigint_to_hex(const BigInt *num, char *hex_string);

// --- ECPoint 转换 ---
void point_to_compressed_hex(const ECPoint *P, char *hex_string);
void point_to_uncompressed_hex(const ECPoint *P, char *hex_string);
#endif // SECP256K1_H
