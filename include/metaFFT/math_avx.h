/*
 * Copyright (C) 2014 Peter Kümmel. All rights reserved.
 *
 * This file is part of metaFFT, distributed under the GNU GPL v2 with
 * a Linking Exception. For full terms see the included COPYING file.
 */
#pragma once

#include "math.h"

#include <immintrin.h>
#include <avxintrin.h>


namespace metaFFT
{

    namespace avx
    {

        struct Splitted {
            double r[4];
            double i[4];
        };

        struct Reg
        {
            __m256d re;
            __m256d im;
        };

        inline Reg MUL(Reg a, Reg b) {
            return {_mm256_sub_pd(_mm256_mul_pd(a.re,b.re),_mm256_mul_pd(a.im,b.im)),
                    _mm256_add_pd(_mm256_mul_pd(a.re,b.im),_mm256_mul_pd(a.im,b.re)) };
        }

        inline Reg MULJ(Reg a, Reg b) {
            return { _mm256_add_pd(_mm256_mul_pd(a.re,b.re),_mm256_mul_pd(a.im,b.im)),
                     _mm256_sub_pd(_mm256_mul_pd(a.im,b.re),_mm256_mul_pd(a.re,b.im)) };
        }

        inline Reg ADD(Reg a, Reg b) {
            return { _mm256_add_pd(a.re,b.re),
                     _mm256_add_pd(a.im,b.im) };
        }

        inline Reg SUB(Reg a, Reg b) {
            return { _mm256_sub_pd(a.re,b.re),
                     _mm256_sub_pd(a.im,b.im) };
        }

        inline Reg ADD_I(Reg a, Reg b) {
            return { _mm256_sub_pd(a.re,b.im),
                     _mm256_add_pd(a.im,b.re) };
        }

        inline Reg SUB_I(Reg a, Reg b) {
            return { _mm256_add_pd(a.re,b.im),
                    _mm256_sub_pd(a.im,b.re) };
        }

        inline Reg LOAD(double* a) {
            return { _mm256_load_pd(a),
                     _mm256_load_pd(a+4) };
        }

        inline void STORE(double* a, Reg r) {
            _mm256_store_pd(a, r.re);
            _mm256_store_pd(a+4, r.im);
        }

        inline void STORE_INTERLEAVED(double *a, Reg r) {
            _mm256_store_pd(a, _mm256_unpacklo_pd(r.re, r.im));
            _mm256_store_pd(a+4, _mm256_unpackhi_pd(r.re, r.im));
        }

    }

}
