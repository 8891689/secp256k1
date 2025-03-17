/* gcc -O3 -o test test.c secp256k1.c */
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include "secp256k1.h"

// 辅助函数：将 BigInt 转换为压缩的公钥（十六进制字符串）
// 注意：这里假设 secp256k1.h 中已经定义了 BigInt 结构体，且其成员为 data
void bigint_to_compressed_pubkey(const BigInt *x, const BigInt *y, char *pubkey_hex) {
    // 确定前缀 (偶数 y 为 0x02, 奇数 y 为 0x03)
    char prefix = (y->data[0] & 1) ? 0x03 : 0x02; // 使用 .data 访问

    // 将 x 坐标转换为十六进制字符串
    char x_hex[65]; // 64 个十六进制字符 + 空终止符
    bigint_to_hex(x, x_hex); // 调用 bigint_to_hex

    // 创建压缩的公钥字符串
    snprintf(pubkey_hex, 67, "%02X%s", prefix, x_hex); // 2 (前缀) + 64 (x) = 66 + 空终止符
}

int main() {
    // 初始化 secp256k1 参数
    BigInt p, Gx, Gy;
    hex_to_bigint("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", &p);
    hex_to_bigint("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", &Gx);
    hex_to_bigint("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", &Gy);

    ECPoint G;
    G.x = Gx;
    G.y = Gy;
    G.infinity = false;

    // 测试：计算私钥 1 到 n 的公钥，并测量时间
    uint32_t num_iterations = 1000; // 更改此值以控制迭代次数
    clock_t start_time = clock();

    for (uint32_t i = 1; i <= num_iterations; i++) {
        BigInt priv;
        init_bigint(&priv, i);
        ECPoint pub;
        scalar_multiply(&priv, &G, &p, &pub);
        //  这里不进行打印，只进行计算以测量速度
    }

    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("进行 %u 次迭代的时间: %.4f 秒\n", num_iterations, elapsed_time);
    printf("每秒计算次数: %.2f\n", num_iterations / elapsed_time);

    // 打印最后 5 个压缩的公钥
    printf("\n最后 5 个压缩的公钥:\n");
    for (uint32_t i = num_iterations - 4; i <= num_iterations; i++) {
        BigInt priv;
        init_bigint(&priv, i);
        char priv_hex[65]; // 64 个十六进制字符 + 空终止符
        bigint_to_hex(&priv, priv_hex); // 将私钥转换为十六进制

        ECPoint pub;
        scalar_multiply(&priv, &G, &p, &pub);

        char compressed_pubkey[67];  // 压缩的公钥: 2 (前缀) + 64 (x 坐标) + 1 (空终止符)
        bigint_to_compressed_pubkey(&pub.x, &pub.y, compressed_pubkey); // 计算压缩公钥

        printf("私钥 (十进制): %u, (十六进制): %s\n", i, priv_hex);
        printf("压缩的公钥: %s\n", compressed_pubkey);
    }

    return 0;
}
