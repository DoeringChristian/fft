/*
 * 
 * MIT License
 * 
 * Copyright (c) 2021 DoeringChristian
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef FFT_H
#define FFT_H

#ifndef FFT_COS
#include <math.h>
#define FFT_COS(_x) cos(_x)
#endif

#ifndef FFT_SIN
#include <math.h>
#define FFT_SIN(_x) sin(_x)
#endif

#ifndef FFT_SIZE
#include <stddef.h>
#define FFT_SIZE size_t
#endif

#ifndef FFT_FLOAT
#define FFT_FLOAT float
#endif

static inline FFT_SIZE log2z(size_t x){
    FFT_SIZE i = 0;
    for(i = sizeof(x) * 8 - 1;((x >> i) & 0x01)!= 0x01 && i != 0;i--);
    return i;
}

/*
 * Reverses bits ok length k in n.
 *
 * @param n: FFT_SIZE from which to reverse bits.
 * @param k: Number of bits to reverse.
 * @return: FFT_SIZE with reversed bits.
 */
static inline FFT_SIZE fft_bitrev(size_t n, size_t bitnum){
    FFT_SIZE ret = 0;
    FFT_SIZE i = 0;
    for(i = 0;i < bitnum;i++)
        ret |= (((n >> (bitnum-i-1)) & 0x01) << i);
    return ret;
}

static inline FFT_FLOAT *vec2_cmult(float *dst, const float *src1, const float *src2){
    FFT_FLOAT src1_real = src1[0];
    FFT_FLOAT src1_imag = src1[1];
    FFT_FLOAT src2_real = src2[0];
    FFT_FLOAT src2_imag = src2[1];

    dst[0] = src1_real *  src2_real - src1_imag * src2_imag;
    dst[1] = src1_real * src2_imag + src1_imag * src2_real;
    return dst;
}

/*
 * fft without complex dependencies, seperate real and imaginary parts.
 *
 * @param dst_real: Destination vector for real values.
 * @param dst_imag: Destination vector for imaginary values.
 * @param src_real: Source vector of real values.
 * @param src_imag: Source vector of imaginary values.
 * @param n: dimension of the vector.
 * @return: dst_real.
 */
FFT_FLOAT *fftri(float *dst_real, float *dst_imag, const FFT_FLOAT *src_real, const FFT_FLOAT *src_imag, FFT_SIZE n){
    FFT_SIZE log2n = log2z(n);
    n = (1 << log2n);


    for(FFT_SIZE i = 0;i < n;i++){
        dst_real[i] = src_real[fft_bitrev(i, log2n)];
        if(src_imag != NULL)
            dst_imag[i] = src_imag[fft_bitrev(i, log2n)];
        else
            dst_imag[i] = 0;
    }

    for(FFT_SIZE i = 0; i < log2n; i++){
        FFT_SIZE m = 1 << i; //m = 2^(recursion_layer)
        FFT_SIZE m2 = m << 1; // m2 = 2^(recursion_layer + 1)
        //FFT_FLOAT complex w = 1;


        //FFT_FLOAT complex wm = cexp(I * (M_PI / m));
        FFT_FLOAT wm[2] = {cos(M_PI / m), sin(M_PI / m)};
        FFT_FLOAT w[2] = {1.0, 0.0};
        for(FFT_SIZE j = 0;j < m;j++){
            for(FFT_SIZE k = j;k < n;k += m2){

                // k + m is the index of the right element.

                //FFT_FLOAT complex t = w * dst[k + m]; // t = e^(-i*pi*section_elemt/(2^recursion_layer)) * x_right
                FFT_FLOAT t[2] = {
                    w[0] * dst_real[k+m] - w[1] * dst_imag[k+m],
                    w[0] * dst_imag[k+m] + w[1] * dst_real[k+m]
                };
                //FFT_FLOAT complex u = dst[k]; // u = x_left
                FFT_FLOAT u[2] = {
                    dst_real[k],
                    dst_imag[k],
                };
                
                //dst[k] = u + t; // x_left,new = u+t
                dst_real[k] = u[0] + t[0];
                dst_imag[k] = u[1] + t[1];

                //dst[k + m] = u - t; // x_right,new = u-t
                dst_real[k+m] = u[0] - t[0];
                dst_imag[k+m] = u[1] - t[1];

            }
            vec2_cmult(w, w, wm);
        }
    }
    return dst_real;
}

/*
 * fft without complex dependencies.
 *
 * @param dst: A Nx2 matrix where the first row denotes the real numbers and the second the imaginary.
 * @param src: A Nx2 matrix where the first row denotes the real numbers and the second the imaginary.
 * @param n: Number of columns
 * @return: dst.
 */
FFT_FLOAT *fftn2(float *dst, const float *src, FFT_SIZE n){
    return fftri(&dst[0], &dst[n], &src[0], &src[n], n);
}

/*
 * fft without complex dependencies.
 *
 * @param dst: A 2xN matrix where every row is made of the real and imaginary part.
 * @param src: A 2xN matrix where every row is made of the real and imaginary part.
 * @param n: number of rows.
 * @return: dst.
 */
FFT_FLOAT *fft2n(float *dst, const float *src, FFT_SIZE n){
    FFT_SIZE log2n = log2z(n);
    n = (1 << log2n);


    for(FFT_SIZE i = 0;i < n;i++){
        dst[i*2 +0] = src[fft_bitrev(i, log2n)*2 +0]; // real part
        dst[i*2 +1] = src[fft_bitrev(i, log2n)*2 +1]; // img part
    }

    for(FFT_SIZE i = 0; i < log2n; i++){
        FFT_SIZE m = 1 << i; //m = 2^(recursion_layer)
        FFT_SIZE m2 = m << 1; // m2 = 2^(recursion_layer + 1)
        //FFT_FLOAT complex w = 1;


        //FFT_FLOAT complex wm = cexp(I * (M_PI / m));
        FFT_FLOAT wm[2] = {cos(M_PI / m), sin(M_PI / m)};
        FFT_FLOAT w[2] = {1.0, 0.0};
        for(FFT_SIZE j = 0;j < m;j++){
            for(FFT_SIZE k = j;k < n;k += m2){

                // k + m is the index of the right element.

                //FFT_FLOAT complex t = w * dst[k + m]; // t = e^(-i*pi*section_elemt/(2^recursion_layer)) * x_right
                FFT_FLOAT t[2];
                vec2_cmult(t, w, &dst[(k+m) * 2]);
                //FFT_FLOAT complex u = dst[k]; // u = x_left
                FFT_FLOAT u[2] = {
                    dst[k * 2 +0],
                    dst[k * 2 +1]
                };
                
                //dst[k] = u + t; // x_left,new = u+t
                dst[k*2 + 0] = u[0] + t[0];
                dst[k*2 + 1] = u[1] + t[1];

                //dst[k + m] = u - t; // x_right,new = u-t
                dst[(k+m) *2 +0] = u[0] - t[0];
                dst[(k+m) *2 +1] = u[1] - t[1];

            }
            vec2_cmult(w, w, wm);
        }
    }
    return dst;
}

#ifndef __STDC_NO_COMPLEX__
#include <complex.h>

/*
 * Calculates the complex fourier series of a vector.
 *
 * This function only calculates the fourier series on the first 2^floor(log2(n)) elements.
 *
 * @param dst: Destination vector to store the result.
 * @param src: Source vector.
 * @param n: Number of elemnts.
 * @return: dst.
 */
FFT_FLOAT complex *fftc(float complex *dst, const float complex *src, FFT_SIZE n){
    FFT_SIZE log2n = log2z(n);
    n = (1 << log2n);


    for(FFT_SIZE i = 0;i < n;i++)
        dst[i] = src[fft_bitrev(i, log2n)];

    for(FFT_SIZE i = 0; i < log2n; i++){
        FFT_SIZE m = 1 << i; //m = 2^(recursion_layer)
        FFT_SIZE m2 = m << 1; // m2 = 2^(recursion_layer + 1)
        FFT_FLOAT complex w = 1;


        FFT_FLOAT complex wm = cexp(I * (M_PI / m));
        for(FFT_SIZE j = 0;j < m;j++){
            for(FFT_SIZE k = j;k < n;k += m2){

                // k + m is the index of the right element.

                FFT_FLOAT complex t = w * dst[k + m]; // t = e^(-i*pi*section_elemt/(2^recursion_layer)) * x_right
                FFT_FLOAT complex u = dst[k]; // u = x_left

                dst[k] = u + t; // x_left,new = u+t
                dst[k + m] = u - t; // x_right,new = u-t
            }
            w *= wm;
        }
    }
    return dst;
}

#endif

#endif //FFT_H
