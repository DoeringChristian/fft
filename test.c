#include "fft.h"
#include <stdio.h>

static inline int veccn_print(const float complex *src, size_t n){
    printf("(");
    for(size_t i = 0;i < n-1;i++)
        printf("%f + %fi, ", creal(src[i]), cimag(src[i]));
    printf("%f + %fi)\n", creal(src[n-1]), cimag(src[n-1]));
    return 1;
}
static inline int vecn_print(const float *src, size_t n){
    printf("(");
    for(size_t i = 0;i < n-1;i++)
        printf("%f, ", src[i]);
    printf("%f)\n", src[n-1]);
    return 1;
}
static inline int matnm_print(const float *src, size_t n, size_t m){
    for(size_t i = 0;i < m;i++){
        printf("|");
        for(size_t j = 0;j < n-1;j++)
            printf("%f, ", src[j + i * n]);
        printf("%f", src[n-1 + i * n]);
        printf("|\n");
    }
    return 1;
}



void test_log2z(){
    printf("%zu\n", log2z(9));
    printf("%zu\n", log2z(10));
    printf("%zu\n", log2z(7));
    printf("%zu\n", log2z(8));
    printf("%zu\n", log2z(17));
}

void test_fft_bitrev(size_t k){
    printf("Test: fft_bitrev(%zu):\n", k);
    for(size_t n = 0;n < pow(2,k);n++){
        printf("%zu\n", fft_bitrev(n, k));
    }
    printf("\n");
}

void test_fftc(){
    printf("Test: fftc:\n");
    float complex src1[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    float complex dst1[8];

    fftc(dst1, src1, 8);

    veccn_print(dst1, 8);
    printf("\n");

    fftc(src1, dst1, 8);

    veccn_print(src1, 8);

    printf("\n");
}

void test_fft2n(){
    printf("Test: fft2n:\n");
    float src1[16] = {
        1, 0,
        2, 0,
        3, 0,
        4, 0,
        5, 0,
        6, 0,
        7, 0,
        8, 0
    };
    float dst1[16];

    fft2n(dst1, src1, 8);

    matnm_print(dst1, 2, 8);
    printf("\n");

    fft2n(src1, dst1, 8);

    matnm_print(src1, 2, 8);

    printf("\n");
}

void test_fftn2(){
    printf("Test: fftn2:\n");
    float src1[16] = {
        1, 2, 3, 4, 5, 6, 7, 8,
        0, 0, 0, 0, 0, 0, 0, 0,
    };
    float dst1[16];

    fftn2(dst1, src1, 8);

    matnm_print(dst1, 8, 2);
    printf("\n");

    fftn2(src1, dst1, 8);

    matnm_print(src1, 8, 2);

    printf("\n");
}

int main(){
#if 0
    test_fft_bitrev(1);
    test_fft_bitrev(2);
    test_fft_bitrev(3);
    test_fft_bitrev(4);
#endif

    test_log2z();

    test_fftc();
    test_fft2n();
    test_fftn2();
    return 0;
}

