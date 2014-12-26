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

#define MM_PERMUTE(f0,f1,f2,f3) (f1 | f1<<1 | f2<<2 | f3<<3)


namespace metaFFT
{
    namespace split_radix
    {

        __inline __m256d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
        _mm256_loadu2_m128d(double const *addr_hi, double const *addr_lo)
        {
          struct __loadu_pd {
            __m128d v;
          } __attribute__((__packed__, __may_alias__));

          __m256d v256 = _mm256_castpd128_pd256(((struct __loadu_pd*)addr_lo)->v);
          return _mm256_insertf128_pd(v256, ((struct __loadu_pd*)addr_hi)->v, 1);
        }

        __inline __m256d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
        _mm256_loadu2_m128d_fill(double const *addr)
        {
          struct __loadu_pd {
            __m128d v;
          } __attribute__((__packed__, __may_alias__));

          __m256d v256 = _mm256_castpd128_pd256(((struct __loadu_pd*)addr)->v);
          return _mm256_insertf128_pd(v256, ((struct __loadu_pd*)addr)->v, 1);
        }



        namespace avx
        {
            namespace unrolled_loop
            {

                template<int Sign, unsigned N, class C>
                struct butterfly_policy
                {
                    typedef typename C::value_type V;

                    template<unsigned K, unsigned End>
                    struct remaining
                    {
                        static void steps(C* out)
                        {

                            V* d = (V*)&out[0];

                            constexpr C cw1 = metaFFT::polar<2*K, N, C, Sign == -1>();
                            constexpr C cw3 = metaFFT::polar<2*3*K, N, C, Sign == -1>();


                            // https://software.intel.com/sites/landingpage/IntrinsicsGuide/

                            //const V w1Zk_re = w1.real()*Zk_re - w1.imag()*Zk_im;
                            //const V w1Zk_im = w1.imag()*Zk_re + w1.real()*Zk_im;
                            constexpr double w1arr[4] = { cw1.real(), -cw1.imag(), cw3.imag(), cw3.real() };
                            __m256d w1 = _mm256_loadu_pd (w1arr);                   // r,-i,i,r
                            __m256d zk = _mm256_loadu2_m128d_fill(&d[2*(K+N/2)]);   // r, i,r,i
                            __m256d w1Zk = _mm256_mul_pd(w1, zk);                   // r*r,-i*i,i*r,r*i

                            //const V w3Zdk_re = w3.real()*Zdk_re - w3.imag()*Zdk_im;
                            //const V w3Zdk_im = w3.imag()*Zdk_re + w3.real()*Zdk_im;
                            constexpr double w3arr[4] = { cw3.real(), -cw3.imag(), cw3.imag(), cw3.real() };
                            __m256d w3 = _mm256_loadu_pd(w3arr);
                            __m256d zdk = _mm256_loadu2_m128d_fill(&d[2*(K+3*N/4)]);
                            __m256d w3Zdk = _mm256_mul_pd(w3, zdk);

                            __m256d w = _mm256_hadd_pd(w1Zk, w3Zdk); // -> w1_re, w3_re, w1_im, w3_im


                            __m256d w13add = _mm256_hadd_pd(w, w); // -> w1Zk_re+w3Zdk_re, w1Zk_re+w3Zdk_re, w1Zk_im+w3Zdk_im, w1Zk_im+w3Zdk_im

                            //d[2*(K)]       = Uk_re  + (w1Zk_re+w3Zdk_re);
                            //d[2*(K)+1]     = Uk_im  + (w1Zk_im+w3Zdk_im);

                            //d[2*(K+N/2)]   = Uk_re  - (w1Zk_re+w3Zdk_re);
                            //d[2*(K+N/2)+1] = Uk_im  - (w1Zk_im+w3Zdk_im);

                            __m256d uk = _mm256_loadu2_m128d_fill(&d[2*(K)]); // -> r,i,r,i
                            uk =  _mm256_shuffle_pd(uk, uk, MM_PERMUTE(0,0,1,1)); // -> r,r,i,i

                            // Uk_re             Uk_re             Uk_im             Uk_im
                            // -                 +                 -                 +
                            // w1Zk_re+w3Zdk_re, w1Zk_re+w3Zdk_re, w1Zk_im+w3Zdk_im, w1Zk_im+w3Zdk_im
                            __m256d res = _mm256_addsub_pd(uk, w13add); // -> Uk_re-re, Uk_re+re, Uk_im-im, Uk_im+im


                            __m256d w13sub = _mm256_hsub_pd(w, w); // -> w1Zk_re-w3Zdk_re, w1Zk_re-w3Zdk_re, w1Zk_im-w3Zdk_im, w1Zk_im-w3Zdk_im

                            //d[2*(K+N/4)]     = Uk2_re - Sign*(+(w1Zk_im-w3Zdk_im));
                            //d[2*(K+N/4)+1]   = Uk2_im + Sign*(+(w1Zk_re-w3Zdk_re));

                            //d[2*(K+3*N/4)]   = Uk2_re + Sign*(+(w1Zk_im-w3Zdk_im));
                            //d[2*(K+3*N/4)+1] = Uk2_im - Sign*(+(w1Zk_re-w3Zdk_re));

                            //const float64x2 Uk2 = load_splat<float64x2>(&d[2*(K+N/4)]);
                            __m256d uk2 = _mm256_loadu2_m128d_fill(&d[2*(K+N/4)]);  // -> r,i,r,i
                            uk2 = _mm256_shuffle_pd(uk2, uk2, MM_PERMUTE(1,1,0,0));      // -> i,i,r,r

                            // Uk2_im            Uk2_im            Uk2_re            Uk2_re
                            // -                 +                 -                 +
                            // w1Zk_re+w3Zdk_re, w1Zk_re+w3Zdk_re, w1Zk_im+w3Zdk_im, w1Zk_im+w3Zdk_im

                            __m256d res2 = _mm256_addsub_pd(uk2, w13sub); // -> Uk2_im-re, Uk2_im+re, Uk2_re-im, Uk2_re+im


                            double duk[4];
                            double duk2[4];

                            _mm256_storeu_pd((double*)duk, res);
                            _mm256_storeu_pd((double*)duk2, res2);

                            d[2*(K)]       = duk[1];
                            d[2*(K)+1]     = duk[3];
                            d[2*(K+N/2)]   = duk[0];
                            d[2*(K+N/2)+1] = duk[2];

                            d[2*(K+N/4)]     = duk2[2];
                            d[2*(K+N/4)+1]   = duk2[1];
                            d[2*(K+3*N/4)]   = duk2[3];
                            d[2*(K+3*N/4)+1] = duk2[0];

                            remaining<K+1, End>::steps(out);

                        }

                    };

                    template<unsigned End>
                    struct remaining<End, End>
                    {
                        static void steps(C*)
                        {}
                    };


                    static void loop(C* out)
                    {
                        remaining<0, N/4>::steps(out);
                    }

                };

            }

        }

    }
}


