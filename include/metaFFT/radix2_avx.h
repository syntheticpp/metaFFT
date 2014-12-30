/*
 * Copyright (C) 2014 Peter Kümmel. All rights reserved.
 *
 * This file is part of metaFFT, distributed under the GNU GPL v2 with
 * a Linking Exception. For full terms see the included COPYING file.
 */
 #pragma once

#include "math_avx.h"



namespace metaFFT
{
    namespace radix2
    {
        namespace avx
        {

            template<unsigned N, class C>
            struct bit_reverse_policy
            {
                typedef typename C::value_type V;

                static void bit_reverse(C* data)
                {
                    V* d = (V*)&data[0];
                    unsigned j = 1;
                    for (unsigned i = 1; i < 2*N; i += 2) {
                        if (j>i) {
                            std::swap(d[j-1], d[i-1]);
                            std::swap(d[j],   d[i]);
                        }
                        unsigned m = N;
                        while (m>=2 && j>m) {
                            j -= m;
                            m >>= 1;
                        }
                        j += m;
                    }
                }
            };


            template<int Sign, unsigned N, class C, bool Unpack>
            struct butterfly_policy
            {
                typedef typename C::value_type V;

                template<unsigned K, unsigned End>
                struct remaining
                {
                    static void steps(C* out)
                    {
                        constexpr C w0 = metaFFT::polar<K+0, N, C, Sign == -1>();
                        constexpr C w1 = metaFFT::polar<K+1, N, C, Sign == -1>();
                        constexpr C w2 = metaFFT::polar<K+3, N, C, Sign == -1>();
                        constexpr C w3 = metaFFT::polar<K+4, N, C, Sign == -1>();
                        constexpr double w[] = { w0.real(), w2.real(), w1.real(), w3.real(),
                                                 w0.imag(), w2.imag(), w1.imag(), w3.imag(), };

                        __m256d w_re = _mm256_load_pd((double*)&w[0]);
                        __m256d w_im = _mm256_load_pd((double*)&w[4]);

                        if (Unpack) {
                            //TODO
                        } else {
                            __m256d Ok_re = _mm256_load_pd((double*)&out[K+N/2]);
                            __m256d Ok_im = _mm256_load_pd((double*)&out[K+N/2+2]);
                            __m256d Ek_re = _mm256_load_pd((double*)&out[K]);
                            __m256d Ek_im = _mm256_load_pd((double*)&out[K+2]);
                            __m256d wOk_re = _mm256_sub_pd(_mm256_mul_pd(Ok_re,w_re),_mm256_mul_pd(Ok_im,w_im));
                            __m256d wOk_im = _mm256_add_pd(_mm256_mul_pd(Ok_re,w_im),_mm256_mul_pd(Ok_im,w_re));
                            _mm256_store_pd((double*)(out+K), _mm256_add_pd(Ek_re, wOk_re));
                            _mm256_store_pd((double*)(out+K+2), _mm256_add_pd(Ek_im, wOk_im));
                            _mm256_store_pd((double*)(out+K+N/2), _mm256_sub_pd(Ek_re, wOk_re));
                            _mm256_store_pd((double*)(out+K+N/2+2), _mm256_sub_pd(Ek_im, wOk_im));
                            static_assert(K+N/2+2 <= N, "wrong K");
                        }

                        remaining<K+2, End>::steps(out);
                    }
                };

                template<unsigned End>
                struct remaining<End, End> { static void steps(C*) {} };


                static void loop(C* cplx)
                {
                    remaining<0, N/2>::steps(cplx);
                }
            };

        }
    }
}
