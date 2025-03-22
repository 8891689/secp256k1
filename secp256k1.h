/*Author: 8891689
  Assist in creation ：ChatGPT 
 */
#ifndef SECP256K1_H
#define SECP256K1_H

#include <stdint.h>
#include <stdbool.h>

#define BIGINT_SIZE 8 // 定义 BigInt 的大小

// 256位大整数（8个32位无符号整数，小端存储：data[0] 为最低位）
typedef struct {
    uint32_t data[8];
} BigInt;

// ECC 点（仿射坐标）
typedef struct {
    BigInt x;
    BigInt y;
    bool infinity;  // true 表示无穷远点
} ECPoint;

// 大整数基础运算声明
void init_bigint(BigInt *x, uint32_t val);
void copy_bigint(BigInt *dest, const BigInt *src);
int compare_bigint(const BigInt *a, const BigInt *b);
bool is_zero(const BigInt *a);
int get_bit(const BigInt *a, int i);

void ptx_u256Add(BigInt *res, const BigInt *a, const BigInt *b);
void ptx_u256Sub(BigInt *res, const BigInt *a, const BigInt *b);

// 9 字（32 位）扩展运算声明（用于模乘中中间结果的折叠）
void multiply_bigint_by_const(const BigInt *a, uint32_t c, uint32_t result[9]);
void shift_left_word(const BigInt *a, uint32_t result[9]);
void add_9word(uint32_t r[9], const uint32_t addend[9]);
void convert_9word_to_bigint(const uint32_t r[9], BigInt *res);

// 模乘及模约简相关函数
void mul_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p);
void efficient_mod(BigInt *r, const BigInt *a, const BigInt *p);
void mod_generic(BigInt *r, const BigInt *a, const BigInt *p);
void sub_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p);
void add_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p);
void modexp(BigInt *res, const BigInt *base, const BigInt *exp, const BigInt *p);
void mod_inverse(BigInt *res, const BigInt *a, const BigInt *p);

// ECC 点运算函数
void point_set_infinity(ECPoint *P);
void point_copy(ECPoint *dest, const ECPoint *src);
void point_add(ECPoint *R, const ECPoint *P, const ECPoint *Q, const BigInt *p);
void double_point(ECPoint *R, const ECPoint *P, const BigInt *p);
void scalar_multiply(const BigInt *d, const ECPoint *G, const BigInt *p, ECPoint *result);

// 辅助工具函数
void print_bigint(const BigInt *b);
void hex_to_bigint(const char* hex, BigInt *b);
void bigint_to_hex(const BigInt *num, char *hex_string); // 函数原型

#endif // SECP256K1_H

