#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include "secp256k1.h"

// ------------------------------
// Supplement missing macro definitions (if secp256k1.h is not defined)
// ------------------------------
#ifndef BIGINT_WORDS
#define BIGINT_WORDS 8  // Adjust according to the actual definition
#endif

// ------------------------------
// Optimized helper function
// ------------------------------
void bigint_to_compressed_pubkey(const BigInt *x, const BigInt *y, char *pubkey_hex) {
    // Directly access the least significant bit to determine parity (data[0] is the least significant bit in little-endian)
    char prefix = (y->data[0] & 1) ? 0x03 : 0x02;

    // Optimize memory access order
    char x_hex[65] = {0};
    for (int i = BIGINT_WORDS-1; i >= 0; i--) {
        sprintf(x_hex + (BIGINT_WORDS-1-i)*8, "%08x", x->data[i]);
    }

    snprintf(pubkey_hex, 67, "%02x%s", prefix, x_hex);
}

// ------------------------------
// Main test function (corrected version)
// ------------------------------
int main() {
    // Initialize secp256k1 parameters (little-endian storage)
    BigInt p, Gx, Gy;
    hex_to_bigint("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", &p);
    hex_to_bigint("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", &Gx);
    hex_to_bigint("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", &Gy);

    // Initialize base point G in Jacobian coordinates (Z=1)
    ECPointJac G_jac;
    init_point_jac(&G_jac, false);
    copy_bigint(&G_jac.X, &Gx);
    copy_bigint(&G_jac.Y, &Gy);
    init_bigint(&G_jac.Z, 1);

    // Warm up cache
    ECPointJac temp_jac;
    BigInt test_priv;
    init_bigint(&test_priv, 1);
    scalar_multiply_jac(&temp_jac, &G_jac, &test_priv, &p);


    // Performance test
    const uint32_t num_iterations = 40000;
    clock_t start = clock();

    for (uint32_t i = 1; i <= num_iterations; ++i) {
        BigInt priv;
        init_bigint(&priv, i);

        ECPointJac result_jac;
        scalar_multiply_jac(&result_jac, &G_jac, &priv, &p);
    }

    double elapsed = (double)(clock() - start)/CLOCKS_PER_SEC;
    printf("[Performance] %u iterations took: %.3fs (%.1f ops/s)\n\n",
           num_iterations, elapsed, num_iterations/elapsed);

    // Verify the last 5 results
    printf("Verification of the last 5 public keys:\n");
    for (uint32_t i = num_iterations-4; i <= num_iterations; ++i) {
        BigInt priv;
        init_bigint(&priv, i);

        ECPointJac result_jac;
        scalar_multiply_jac(&result_jac, &G_jac, &priv, &p);

        ECPoint pub;
        jacobian_to_affine(&pub, &result_jac, &p);

        char priv_hex[65] = {0}, pub_hex[67] = {0};
        bigint_to_hex(&priv, priv_hex);
        bigint_to_compressed_pubkey(&pub.x, &pub.y, pub_hex);

        printf("Private key %u:\n  Hex: %s\n  Public Key: %s\n", i, priv_hex, pub_hex);
    }

    return 0;
}
