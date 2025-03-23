/* gcc -O3 -o test test.c secp256k1.c */
#include "secp256k1.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

int main() {
    // Initialize secp256k1 parameters
    const char *p_hex = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";
    BigInt p;
    hex_to_bigint(p_hex, &p);

    // Initialize base point G
    const char *Gx_hex = "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798";
    const char *Gy_hex = "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8";
    ECPoint G;
    hex_to_bigint(Gx_hex, &G.x);
    hex_to_bigint(Gy_hex, &G.y);
    G.infinity = false;

    // Convert to Jacobian coordinates
    ECPointJac G_jac;
    copy_bigint(&G_jac.X, &G.x);
    copy_bigint(&G_jac.Y, &G.y);
    init_bigint(&G_jac.Z, 1);
    G_jac.infinity = false;

    // Store last 5 public keys
    ECPoint last_5[5];
    int count = 20000;
    int start_index = count - 5;

    clock_t start = clock();

    for(int k = 1; k <= count; k++) {
        BigInt scalar;
        init_bigint(&scalar, k);

        ECPointJac result_jac;
        scalar_multiply_jac(&result_jac, &G_jac, &scalar, &p);

        ECPoint result_affine;
        jacobian_to_affine(&result_affine, &result_jac, &p);

        // Store last 5 keys
        if(k > start_index) {
            int idx = k - start_index - 1;
            copy_bigint(&last_5[idx].x, &result_affine.x);
            copy_bigint(&last_5[idx].y, &result_affine.y);
            last_5[idx].infinity = result_affine.infinity;
        }

        // Free memory
        free_bigint(&scalar);
        free_bigint(&result_jac.X);
        free_bigint(&result_jac.Y);
        free_bigint(&result_jac.Z);
        free_bigint(&result_affine.x);
        free_bigint(&result_affine.y);
    }

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    printf("Generated %d public keys in %.4f seconds\n", count, time_spent);
    printf("Last 5 public keys and corresponding private keys:\n");

    char compressed_hex[66] = {0};
    char uncompressed_hex[130] = {0};
    char privkey_hex[65] = {0};

    for(int i = 0; i < 5; i++) {
        int current_k = start_index + i + 1;

        // Convert private key to hex
        BigInt privkey_bigint;
        init_bigint(&privkey_bigint, current_k);
        bigint_to_hex(&privkey_bigint, privkey_hex);
        free_bigint(&privkey_bigint);

        // Generate public key formats
        point_to_compressed_hex(&last_5[i], compressed_hex);
        point_to_uncompressed_hex(&last_5[i], uncompressed_hex);

        // Print results
        printf("Private Key %d: %s\n", current_k, privkey_hex);
        printf("Compressed Public Key: %s\n", compressed_hex);
        printf("Uncompressed Public Key: %s\n\n", uncompressed_hex);
    }

    // Cleanup
    free_bigint(&p);
    free_bigint(&G.x);
    free_bigint(&G.y);
    free_bigint(&G_jac.X);
    free_bigint(&G_jac.Y);
    free_bigint(&G_jac.Z);
    for(int i = 0; i < 5; i++) {
        free_bigint(&last_5[i].x);
        free_bigint(&last_5[i].y);
    }

    return 0;
}

