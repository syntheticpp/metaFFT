/*
 * Copyright (C) 2014 Peter Kümmel. All rights reserved.
 *
 * This file is part of metaFFT, distributed under the GNU GPL v2 with
 * a Linking Exception. For full terms see the included COPYING file.
 */
 #pragma once

#include "math.h"

#include <immintrin.h>


namespace metaFFT
{
    namespace radix2
    {
        namespace ctran
        {
            namespace sse2
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
                            static void steps(V* d)
                            {
                                constexpr C b = metaFFT::polar<K, N, C, Sign == -1>();

                                __m128d ta, ta2, tb, rs; // meaning of chars: t=tmp, rs=result, r/i=real/imag

                                // mul
                                
                                // (ar + Iai) * (br + Ibi) = (1.)ar*br + I(2.)ar*bi + I(3.)ai*br + I*I(4.)ai*bi
                                //                         = (1.)ar*br -  (4.)ai*bi + I(2.)ar*bi +   I(3.)ai*br
                                ta = _mm_loaddup_pd(&d[K+N]);        // -> {ar, ar}
                                ta2 = _mm_loaddup_pd(&d[K+N+1]);     // -> {ai, ai}
                                tb = _mm_set_pd(b.imag(), b.real()); // -> {br, bi}
                                rs = _mm_mul_pd(ta, tb);             // -> {(1.)ar*br, (2.)ar*bi}

                                tb = _mm_shuffle_pd(tb, tb, 1);      // -> {bi, br}
                                ta = _mm_mul_pd(tb, ta2);             // -> {(4.)ai*bi, (3.)ai*br}

                                rs = _mm_addsub_pd(rs, ta2);          // ->  {(1.) - (4.), (2.) + (3.)}

                                
                                // add/sub

                                //d[K+N]   = d[K]   - t.real;
                                //d[K+N+1] = d[K+1] - t.imag;
                                ta = _mm_set_pd(d[K+1], d[K]);
                                tb = _mm_sub_pd(ta, rs);

                                //d[K]     += t.real;
                                //d[K+1]   += t.imag;
                                ta = _mm_add_pd(ta, rs);
                                _mm_storeu_pd(&d[K+N], tb);
                                _mm_storeu_pd(&d[K], ta);

                                remaining<K+2, End>::steps(d); 

                            }
                        };

                        template<unsigned End>
                        struct remaining<End, End> { static void steps(V*) {} };


                        static void loop(C* data)
                        {
                            V* d = (V*)&data[0];
                            remaining<0, N>::steps(d);
                        }
                    };
                }
            }
        }
    }
}
