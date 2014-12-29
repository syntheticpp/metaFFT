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
    namespace split_radix
    {

        namespace avx
        {

            namespace unrolled_loop
            {

                template<int Sign, unsigned N, class C, bool Unpack>
                struct butterfly_policy
                {
                    typedef typename C::value_type V;


                    template<unsigned K, unsigned End>
                    struct remaining
                    {

                        static void steps(C* out)
                        {
                            using namespace metaFFT::avx;

                            constexpr C w1_0 = metaFFT::polar<2*(K+0), N, C, Sign == -1>();
                            constexpr C w1_1 = metaFFT::polar<2*(K+1), N, C, Sign == -1>();
                            constexpr C w1_2 = metaFFT::polar<2*(K+2), N, C, Sign == -1>();
                            constexpr C w1_3 = metaFFT::polar<2*(K+3), N, C, Sign == -1>();

                            constexpr C w3_0 = metaFFT::polar<2*3*(K+0), N, C, Sign == -1>();
                            constexpr C w3_1 = metaFFT::polar<2*3*(K+1), N, C, Sign == -1>();
                            constexpr C w3_2 = metaFFT::polar<2*3*(K+2), N, C, Sign == -1>();
                            constexpr C w3_3 = metaFFT::polar<2*3*(K+3), N, C, Sign == -1>();

                            constexpr Splitted arr[] = {
                                                            {{ w1_0.real(), w1_2.real(), w1_1.real(), w1_3.real() },
                                                             { w1_0.imag(), w1_2.imag(), w1_1.imag(), w1_3.imag() }},
                                                            {{ w3_0.real(), w3_2.real(), w3_1.real(), w3_3.real() },
                                                             { w3_0.imag(), w3_2.imag(), w3_1.imag(), w3_3.imag() }}
                                                           };


                           // always in split format
                            Reg Uk = LOAD((double*)&out[K]);
                            Reg Zk = LOAD((double*)&out[K+N/2]);
                            Reg Uk2 = LOAD((double*)&out[K+N/4]);
                            Reg Zdk = LOAD((double*)&out[K+3*N/4]);
                            Reg w1 = LOAD((double*)&arr[0]);
                            Reg w3 = LOAD((double*)&arr[1]);

                            Reg w3Zdk = MUL(w3, Zdk);
                            Reg w1Zk = MUL(w1, Zk);
                            Reg sum = ADD(w1Zk, w3Zdk);
                            Reg dif = SUB(w1Zk, w3Zdk);


                            if (Unpack) {
                                // switch back to interleaved format
                                STORE_INTERLEAVED((double*)&out[K], ADD(Uk, sum));
                                STORE_INTERLEAVED((double*)&out[K+N/2], SUB(Uk, sum));
                                STORE_INTERLEAVED((double*)&out[K+N/4], SUB_I(Uk2, dif));
                                STORE_INTERLEAVED((double*)&out[K+3*N/4], ADD_I(Uk2, dif));
                            } else {
                                STORE((double*)&out[K], ADD(Uk, sum));
                                STORE((double*)&out[K+N/2], SUB(Uk, sum));
                                STORE((double*)&out[K+N/4], SUB_I(Uk2, dif));
                                STORE((double*)&out[K+3*N/4], ADD_I(Uk2, dif));
                            }

                            remaining<K+4, End>::steps(out);
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


