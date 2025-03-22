/* gcc -O3 -o test test.c secp256k1.c */
#include "secp256k1.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

int main() {
    // 初始化secp256k1参数
    const char *p_hex = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";
    BigInt p;
    hex_to_bigint(p_hex, &p);

    // 初始化基点G
    const char *Gx_hex = "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798";
    const char *Gy_hex = "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8";
    ECPoint G;
    hex_to_bigint(Gx_hex, &G.x);
    hex_to_bigint(Gy_hex, &G.y);
    G.infinity = false;

    // 转换为雅可比坐标
    ECPointJac G_jac;
    copy_bigint(&G_jac.X, &G.x);
    copy_bigint(&G_jac.Y, &G.y);
    init_bigint(&G_jac.Z, 1);
    G_jac.infinity = false;

    // 存储最后5个公钥
    ECPoint last_5[5];
    int count = 50000;
    int start_index = count -5;

    clock_t start = clock();

    for(int k = 1; k <= count; k++) {
        BigInt scalar;
        init_bigint(&scalar, k); // 假设init_bigint可以设置大整数为k的值

        ECPointJac result_jac;
        scalar_multiply_jac(&result_jac, &G_jac, &scalar, &p); // 标量乘法

        ECPoint result_affine;
        jacobian_to_affine(&result_affine, &result_jac, &p);

        // 存储最后5个
        if(k > start_index) {
            int idx = k - start_index -1;
            copy_bigint(&last_5[idx].x, &result_affine.x);
            copy_bigint(&last_5[idx].y, &result_affine.y);
            last_5[idx].infinity = result_affine.infinity;
        }

        // 释放内存，假设有free_bigint函数
        free_bigint(&scalar);
        free_bigint(&result_jac.X);
        free_bigint(&result_jac.Y);
        free_bigint(&result_jac.Z);
        free_bigint(&result_affine.x);
        free_bigint(&result_affine.y);
    }

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    printf("生成%d个公钥耗时: %.4f秒\n", count, time_spent);
    printf("最后5个公钥坐标:\n");

    char hexStr[256] = {0};
    for(int i=0; i<5; i++) {
        memset(hexStr, 0, sizeof(hexStr));
        bigint_to_hex(&last_5[i].x, hexStr);
        printf("公钥%d x: %s\n", count -5 + i +1, hexStr);

        memset(hexStr, 0, sizeof(hexStr));
        bigint_to_hex(&last_5[i].y, hexStr);
        printf("公钥%d y: %s\n", count -5 + i +1, hexStr);
    }

    // 释放其他资源_b_bigint(&p);
    free_bigint(&G.x);
    free_bigint(&G.y);
    free_bigint(&G_jac.X);
    free_bigint(&G_jac.Y);
    free_bigint(&G_jac.Z);
    for(int i=0; i<5; i++) {
        free_bigint(&last_5[i].x);
        free_bigint(&last_5[i].y);
    }

    return 0;
}
