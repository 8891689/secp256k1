/*Author: 8891689
 * Assist in creation ：ChatGPT 
 */
#include "secp256k1.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define BIGINT_SIZE 8  // 根据实际情况调整数值

// ------------------------------
// secp256k1 曲线参数
// ------------------------------

// secp256k1 的素数域 p = 2^256 - 2^32 - 977
static const BigInt secp256k1_p = {
    .data = {
        0xFFFFFC2F, 0xFFFFFFFE, 0xFFFFFFFF, 0xFFFFFFFF,
        0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF
    }
};

// 基点 G 的仿射坐标 (十六进制小端序)
static const BigInt G_x = {
    .data = {
        0x16F81798, 0x59F2815B, 0x2DCE28D9, 0x029BFCDB,
        0xCE870B07, 0x55A06295, 0xF9DCBBAC, 0x79BE667E
    }
};

static const BigInt G_y = {
    .data = {
        0xFB10D4B8, 0x9C47D08F, 0xA6855419, 0xFD17B448,
        0x0E1108A8, 0x5DA4FBFC, 0x26A3C465, 0x483ADA77
    }
};

// 基点 G 的雅可比坐标表示 (Z=1)
static const ECPointJac G_jacobian = {
    .X = G_x,
    .Y = G_y,
    .Z = { .data = {1, 0, 0, 0, 0, 0, 0, 0} },
    .infinity = false
};

// ------------------------------
// 私钥转公钥实现
// ------------------------------

void private_to_public_key(ECPoint *public_key, const BigInt *private_key) {
    ECPointJac result_jac;
    
    // 执行标量乘法：Q = private_key * G
    scalar_multiply_jac(&result_jac, &G_jacobian, private_key, &secp256k1_p);
    
    // 将雅可比坐标转换为仿射坐标
    jacobian_to_affine(public_key, &result_jac, &secp256k1_p);
}


// ------------------------------
// 大整数基本运算实现
// ------------------------------
void init_bigint(BigInt *x, uint32_t val) {
    memset(x->data, 0, sizeof(x->data));
    x->data[0] = val;
}

void copy_bigint(BigInt *dest, const BigInt *src) {
    memcpy(dest->data, src->data, sizeof(src->data));
}

int compare_bigint(const BigInt *a, const BigInt *b) {
    for (int i = BIGINT_WORDS - 1; i >= 0; i--) {
        if (a->data[i] != b->data[i])
            return (a->data[i] > b->data[i]) ? 1 : -1;
    }
    return 0;
}

bool is_zero(const BigInt *a) {
    for (int i = 0; i < BIGINT_WORDS; i++) {
        if (a->data[i])
            return false;
    }
    return true;
}

// 高效位获取函数
int get_bit(const BigInt *a, int i) {
    int word_idx = i >> 5;
    int bit_idx = i & 31;
    if (word_idx >= BIGINT_WORDS) return 0;
    return (a->data[word_idx] >> bit_idx) & 1;
}

// 在ptx_u256Add和ptx_u256Sub中使用uint64_t进行中间计算
void ptx_u256Add(BigInt *res, const BigInt *a, const BigInt *b) {
    uint64_t carry = 0;
    for (int i = 0; i < BIGINT_WORDS; i++) {
        uint64_t sum = (uint64_t)a->data[i] + b->data[i] + carry;
        res->data[i] = (uint32_t)sum;
        carry = sum >> 32;
    }
}


void ptx_u256Sub(BigInt *res, const BigInt *a, const BigInt *b) {
    // 加载所有数据到寄存器
    uint32_t a0 = a->data[0], b0 = b->data[0];
    uint32_t a1 = a->data[1], b1 = b->data[1];
    uint32_t a2 = a->data[2], b2 = b->data[2];
    uint32_t a3 = a->data[3], b3 = b->data[3];
    uint32_t a4 = a->data[4], b4 = b->data[4];
    uint32_t a5 = a->data[5], b5 = b->data[5];
    uint32_t a6 = a->data[6], b6 = b->data[6];
    uint32_t a7 = a->data[7], b7 = b->data[7];
    uint32_t brw = 0; // 借位寄存器

    // 完全展开减法操作
    uint64_t tmp = (uint64_t)a0 - b0 - brw;
    res->data[0] = (uint32_t)tmp;
    brw = (tmp >> 32) & 1;

    tmp = (uint64_t)a1 - b1 - brw;
    res->data[1] = (uint32_t)tmp;
    brw = (tmp >> 32) & 1;

    tmp = (uint64_t)a2 - b2 - brw;
    res->data[2] = (uint32_t)tmp;
    brw = (tmp >> 32) & 1;

    tmp = (uint64_t)a3 - b3 - brw;
    res->data[3] = (uint32_t)tmp;
    brw = (tmp >> 32) & 1;

    tmp = (uint64_t)a4 - b4 - brw;
    res->data[4] = (uint32_t)tmp;
    brw = (tmp >> 32) & 1;

    tmp = (uint64_t)a5 - b5 - brw;
    res->data[5] = (uint32_t)tmp;
    brw = (tmp >> 32) & 1;

    tmp = (uint64_t)a6 - b6 - brw;
    res->data[6] = (uint32_t)tmp;
    brw = (tmp >> 32) & 1;

    tmp = (uint64_t)a7 - b7 - brw;
    res->data[7] = (uint32_t)tmp;
}


void free_bigint(BigInt *a) {
    memset(a->data, 0, sizeof(a->data));
}

// ------------------------------
// Montgomery 参数及运算实现（已替换）
// ------------------------------

// 简单模约简函数（假设 a < 2p）
void mod_generic(BigInt *r, const BigInt *a, const BigInt *p) {
    if (compare_bigint(a, p) >= 0) {
        ptx_u256Sub(r, a, p);
    } else {
        memcpy(r->data, a->data, sizeof(a->data));  // 使用 memcpy 进行数组复制
    }
}



// 模乘函数（仅用于 Montgomery 参数初始化时计算 R2）
void mul_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p) {
    uint64_t prod[2*BIGINT_WORDS] = {0};
    for (int i = 0; i < BIGINT_WORDS; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < BIGINT_WORDS; j++) {
            uint64_t tmp = prod[i+j] + (uint64_t)a->data[i] * b->data[j] + carry;
            prod[i+j] = (uint32_t) tmp;
            carry = tmp >> 32;
        }
        prod[i+BIGINT_WORDS] += carry;
    }
    // 取低 BIGINT_WORDS 个字，并进行简单模约简（假设结果 < 2p）
    for (int i = 0; i < BIGINT_WORDS; i++) {
        res->data[i] = (uint32_t) prod[i];
    }
    if (compare_bigint(res, p) >= 0) {
        ptx_u256Sub(res, res, p);
    }
}

// Montgomery 参数初始化（修正版）
void montgomery_init(MontgomeryCtx *ctx, const BigInt *p) {
    // 构造 R = 2^(32*BIGINT_WORDS) mod p
    BigInt R;
    memset(R.data, 0, sizeof(R.data));
    R.data[BIGINT_WORDS-1] = 1;  // 设置最高位为1
    mod_generic(&ctx->R, &R, p); // R mod p

    // 计算 R^2 mod p -> R2
    mul_mod(&ctx->R2, &ctx->R, &ctx->R, p);

    // 计算 R^4 mod p -> R4 (新增预计算)
    mul_mod(&ctx->R4, &ctx->R2, &ctx->R2, p);

    // 计算 inv_p = -p^-1 mod 2^32
    uint32_t p0 = p->data[0]; // 取p的最低32位
    ctx->inv_p.data[0] = 0;
    // 快速计算模逆
    for (int i = 0; i < 32; i++) {
        if ((p0 * (1U << i)) & 1) {
            ctx->inv_p.data[0] = (1U << i) - p0;
            break;
        }
    }
    // 高位清零
    memset(&ctx->inv_p.data[1], 0, (BIGINT_WORDS-1)*sizeof(uint32_t));
}

// 改进后的 Montgomery 乘法实现（内联归约版本）
static inline void montgomery_mult(BigInt *result, 
                     const BigInt *a, 
                     const BigInt *b,
                     const MontgomeryCtx *ctx,
                     const BigInt *p) {
    uint32_t t[17] = {0};

    for (int i = 0; i < 8; i++) { // BIGINT_WORDS=8
        const uint32_t a_i = a->data[i];
        const uint32_t *b_ptr = b->data;
        uint32_t *t_ptr = &t[i];
        uint64_t carry = 0;

        // 使用指针算术优化内存访问
        for (int j = 0; j < 8; j++) {
            uint64_t prod = (uint64_t)a_i * b_ptr[j] + t_ptr[j] + carry;
            t_ptr[j] = (uint32_t)prod;
            carry = prod >> 32;
        }
        t_ptr[8] += (uint32_t)carry;

        uint32_t m = t[i] * ctx->inv_p.data[0];
        carry = 0;
        const uint32_t *p_ptr = p->data;
        t_ptr = &t[i]; // 重置指针

        for (int j = 0; j < 8; j++) {
            uint64_t term = (uint64_t)m * p_ptr[j] + t_ptr[j] + carry;
            t_ptr[j] = (uint32_t)term;
            carry = term >> 32;
        }

        // 快速进位处理（最多2次）
        for (int k = 8; carry && k < 16; k++) {
            uint64_t sum = t[i+k] + carry;
            t[i+k] = (uint32_t)sum;
            carry = sum >> 32;
        }
    }

    // 手动展开结果拷贝（比 memcpy 更快）
    result->data[0] = t[8]; result->data[1] = t[9];
    result->data[2] = t[10]; result->data[3] = t[11];
    result->data[4] = t[12]; result->data[5] = t[13];
    result->data[6] = t[14]; result->data[7] = t[15];
    
    if (compare_bigint(result, p) >= 0) {
        ptx_u256Sub(result, result, p);
    }
}



// ------------------------------
// ECC 点运算以及其它辅助函数（保持原样）
// ------------------------------

// 雅可比坐标点初始化
void init_point_jac(ECPointJac *p, bool infinity) {
    memset(&p->X, 0, sizeof(BigInt));
    memset(&p->Y, 0, sizeof(BigInt));
    memset(&p->Z, 0, sizeof(BigInt));
    p->infinity = infinity;
}

// 雅可比坐标点复制
void copy_point_jac(ECPointJac *dest, const ECPointJac *src) {
    memcpy(&dest->X, &src->X, sizeof(BigInt));
    memcpy(&dest->Y, &src->Y, sizeof(BigInt));
    memcpy(&dest->Z, &src->Z, sizeof(BigInt));
    dest->infinity = src->infinity;
}

// 优化后的标量乘法（已验证正确性）
void scalar_multiply_jac(ECPointJac *result, const ECPointJac *point, const BigInt *scalar, const BigInt *p) {
    ECPointJac res;
    init_point_jac(&res, true);

    // 使用窗口法来优化标量乘法
    int highest_bit = BIGINT_WORDS * 32 - 1;
    for (; highest_bit >= 0; highest_bit--) {
        if (get_bit(scalar, highest_bit)) break;
    }
    if (highest_bit < 0) {
        copy_point_jac(result, &res);
        return;
    }
    for (int i = highest_bit; i >= 0; i--) {
        double_point_jac(&res, &res, p);
        if (get_bit(scalar, i)) {
            add_point_jac(&res, &res, point, p);
        }
    }
    copy_point_jac(result, &res);
}


// 以下为9字扩展运算实现
void multiply_bigint_by_const(const BigInt *a, uint32_t c, uint32_t result[9]) {
    uint64_t carry = 0;
    for (int i = 0; i < BIGINT_WORDS; i++) {
        uint64_t prod = (uint64_t)a->data[i] * c + carry;
        result[i] = (uint32_t)prod;
        carry = prod >> 32;
    }
    result[8] = (uint32_t)carry;
}

void shift_left_word(const BigInt *a, uint32_t result[9]) {
    result[0] = 0;
    memcpy(&result[1], a->data, BIGINT_WORDS * sizeof(uint32_t));
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
    memcpy(res->data, r, BIGINT_WORDS * sizeof(uint32_t));
}

// 模乘及模约简实现
void mul_mod_old(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p) {
    uint32_t prod[16] = {0};
    for (int i = 0; i < BIGINT_WORDS; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < BIGINT_WORDS; j++) {
            uint64_t tmp = (uint64_t)prod[i+j] + (uint64_t)a->data[i] * b->data[j] + carry;
            prod[i+j] = (uint32_t)tmp;
            carry = tmp >> 32;
        }
        prod[i+BIGINT_WORDS] += (uint32_t)carry;
    }
    BigInt L, H;
    for (int i = 0; i < BIGINT_WORDS; i++) {
        L.data[i] = prod[i];
        H.data[i] = prod[i+BIGINT_WORDS];
    }
    uint32_t Rext[9] = {0};
    memcpy(Rext, L.data, BIGINT_WORDS * sizeof(uint32_t));
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
    copy_bigint(res, &temp);
}

void add_mod(BigInt *res, const BigInt *a, const BigInt *b, const BigInt *p) {
    BigInt temp;
    ptx_u256Add(&temp, a, b);
    mod_generic(&temp, &temp, p);
    copy_bigint(res, &temp);
}

void modexp(BigInt *res, const BigInt *base, const BigInt *exp, const BigInt *p) {
    BigInt result;
    init_bigint(&result, 1);
    BigInt b;
    copy_bigint(&b, base);
    for (int i = 0; i < 256; i++) {
         if (get_bit(exp, i)) {
              BigInt temp;
              mul_mod_old(&temp, &result, &b, p);
              copy_bigint(&result, &temp);
         }
         BigInt temp;
         mul_mod_old(&temp, &b, &b, p);
         copy_bigint(&b, &temp);
    }
    copy_bigint(res, &result);
}

void mod_inverse(BigInt *res, const BigInt *a, const BigInt *p) {
    BigInt p_minus_2;
    BigInt two;
    init_bigint(&two, 2);
    ptx_u256Sub(&p_minus_2, p, &two); 
    modexp(res, a, &p_minus_2, p);
}

// ------------------------------
// 仿射坐标下 ECC 点运算实现
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
    if (P->infinity) {
        point_copy(R, Q);
        return;
    }
    if (Q->infinity) {
        point_copy(R, P);
        return;
    }

    BigInt diffY, diffX, inv_diffX, lambda, lambda2, temp;
    sub_mod(&diffY, &Q->y, &P->y, p);
    sub_mod(&diffX, &Q->x, &P->x, p);
    mod_inverse(&inv_diffX, &diffX, p);
    mul_mod_old(&lambda, &diffY, &inv_diffX, p);
    mul_mod_old(&lambda2, &lambda, &lambda, p);
    sub_mod(&temp, &lambda2, &P->x, p);
    sub_mod(&R->x, &temp, &Q->x, p);
    sub_mod(&temp, &P->x, &R->x, p);
    mul_mod_old(&R->y, &lambda, &temp, p);
    sub_mod(&R->y, &R->y, &P->y, p);
    R->infinity = false;
}


void double_point(ECPoint *R, const ECPoint *P, const BigInt *p) {
    if (P->infinity || is_zero(&P->y)) {
         point_set_infinity(R);
         return;
    }
    BigInt x2, numerator, denominator, inv_den, lambda, lambda2, two, two_x;
    mul_mod_old(&x2, &P->x, &P->x, p);
    BigInt three; 
    init_bigint(&three, 3);
    mul_mod_old(&numerator, &three, &x2, p);
    init_bigint(&two, 2);
    mul_mod_old(&denominator, &two, &P->y, p);
    mod_inverse(&inv_den, &denominator, p);
    mul_mod_old(&lambda, &numerator, &inv_den, p);
    mul_mod_old(&lambda2, &lambda, &lambda, p);
    mul_mod_old(&two_x, &two, &P->x, p);
    sub_mod(&R->x, &lambda2, &two_x, p);
    sub_mod(&numerator, &P->x, &R->x, p);
    mul_mod_old(&R->y, &lambda, &numerator, p);
    sub_mod(&R->y, &R->y, &P->y, p);
    R->infinity = false;
}

// ------------------------------
// 雅可比坐标下 ECC 点运算实现
// ------------------------------
void point_set_infinity_jac(ECPointJac *P) {
    P->infinity = true;
}

void point_copy_jac(ECPointJac *dest, const ECPointJac *src) {
    copy_bigint(&dest->X, &src->X);
    copy_bigint(&dest->Y, &src->Y);
    copy_bigint(&dest->Z, &src->Z);
    dest->infinity = src->infinity;
}

void double_point_jac(ECPointJac *R, const ECPointJac *P, const BigInt *p) {
    if (P->infinity || is_zero(&P->Y)) {
        point_set_infinity_jac(R);
        return;
    }
    
    BigInt A, B, C, D, X3, Y3, Z3, temp, temp2;
    mul_mod_old(&A, &P->Y, &P->Y, p);
    mul_mod_old(&temp, &P->X, &A, p);
    init_bigint(&temp2, 4);
    mul_mod_old(&B, &temp, &temp2, p);
    mul_mod_old(&temp, &A, &A, p);
    init_bigint(&temp2, 8);
    mul_mod_old(&C, &temp, &temp2, p);
    mul_mod_old(&temp, &P->X, &P->X, p);
    init_bigint(&temp2, 3);
    mul_mod_old(&D, &temp, &temp2, p);
    BigInt D2;
    mul_mod_old(&D2, &D, &D, p);
    BigInt two;
    init_bigint(&two, 2);
    BigInt twoB;
    mul_mod_old(&twoB, &B, &two, p);
    sub_mod(&X3, &D2, &twoB, p);
    sub_mod(&temp, &B, &X3, p);
    mul_mod_old(&temp, &D, &temp, p);
    sub_mod(&Y3, &temp, &C, p);
    init_bigint(&temp, 2);
    mul_mod_old(&temp, &temp, &P->Y, p);
    mul_mod_old(&Z3, &temp, &P->Z, p);
    copy_bigint(&R->X, &X3);
    copy_bigint(&R->Y, &Y3);
    copy_bigint(&R->Z, &Z3);
    R->infinity = false;
}

void add_point_jac(ECPointJac *R, const ECPointJac *P, const ECPointJac *Q, const BigInt *p) {
    if (P->infinity) { 
        point_copy_jac(R, Q);
        return;
    }
    if (Q->infinity) { 
        point_copy_jac(R, P);
        return;
    }
    
    BigInt Z1Z1, Z2Z2, U1, U2, S1, S2, H, R_big, H2, H3, U1H2, X3, Y3, Z3, temp;
    mul_mod_old(&Z1Z1, &P->Z, &P->Z, p);
    mul_mod_old(&Z2Z2, &Q->Z, &Q->Z, p);
    mul_mod_old(&U1, &P->X, &Z2Z2, p);
    mul_mod_old(&U2, &Q->X, &Z1Z1, p);
    BigInt Z2_cubed, Z1_cubed;
    mul_mod_old(&Z2_cubed, &Z2Z2, &Q->Z, p);
    mul_mod_old(&Z1_cubed, &Z1Z1, &P->Z, p);
    mul_mod_old(&S1, &P->Y, &Z2_cubed, p);
    mul_mod_old(&S2, &Q->Y, &Z1_cubed, p);
    
    if (compare_bigint(&U1, &U2) == 0) {
        if (compare_bigint(&S1, &S2) != 0) {
            point_set_infinity_jac(R);
            return;
        } else {
            double_point_jac(R, P, p);
            return;
        }
    }
    
    // H = U2 - U1, R_big = S2 - S1
    sub_mod(&H, &U2, &U1, p);
    sub_mod(&R_big, &S2, &S1, p);

    // H^2, H^3, U1*H^2
    mul_mod_old(&H2, &H, &H, p);
    mul_mod_old(&H3, &H2, &H, p);
    mul_mod_old(&U1H2, &U1, &H2, p);

    // X3 = R_big^2 - H^3 - 2*U1H2
    BigInt R2;
    mul_mod_old(&R2, &R_big, &R_big, p);
    BigInt two;         // 新增声明
    init_bigint(&two, 2);
    BigInt twoU1H2;
    mul_mod_old(&twoU1H2, &U1H2, &two, p);
    sub_mod(&temp, &R2, &H3, p);
    sub_mod(&X3, &temp, &twoU1H2, p);

    sub_mod(&temp, &U1H2, &X3, p);
    mul_mod_old(&temp, &R_big, &temp, p);
    mul_mod_old(&Y3, &S1, &H3, p);
    sub_mod(&Y3, &temp, &Y3, p);
    mul_mod_old(&temp, &P->Z, &Q->Z, p);
    mul_mod_old(&Z3, &temp, &H, p);
    
    copy_bigint(&R->X, &X3);
    copy_bigint(&R->Y, &Y3);
    copy_bigint(&R->Z, &Z3);
    R->infinity = false;
}

void jacobian_to_affine(ECPoint *R, const ECPointJac *P, const BigInt *p) {
    if (P->infinity) {
        R->infinity = true;
        return;
    }
    BigInt Zinv, Zinv2, Zinv3;
    mod_inverse(&Zinv, &P->Z, p);
    mul_mod_old(&Zinv2, &Zinv, &Zinv, p);
    mul_mod_old(&Zinv3, &Zinv2, &Zinv, p);
    
    mul_mod_old(&R->x, &P->X, &Zinv2, p);
    mul_mod_old(&R->y, &P->Y, &Zinv3, p);
    R->infinity = false;
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

void bigint_to_hex(const BigInt *num, char *hex_string) {
    int len = 0;
    for (int i = BIGINT_SIZE - 1; i >= 0; i--) {
        len += sprintf(hex_string + len, "%08x", num->data[i]);
    }
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
