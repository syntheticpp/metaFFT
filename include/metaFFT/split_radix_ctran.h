/*
 * Copyright (C) 2014 Peter Kümmel. All rights reserved.
 *
 * This file is part of metaFFT, distributed under the GNU GPL v2 with
 * a Linking Exception. For full terms see the included COPYING file.
 */
#pragma once

#include "math.h"


namespace metaFFT
{
    namespace split_radix
    {
        namespace ctran
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

                            constexpr C w1 = metaFFT::polar<2*K, N, C, Sign == -1>();
                            constexpr C w3 = metaFFT::polar<2*3*K, N, C, Sign == -1>();

                            const V Uk_re  = d[2*(K)];
                            const V Uk_im  = d[2*(K)+1];
                            const V Zk_re  = d[2*(K+N/2)];
                            const V Zk_im  = d[2*(K+N/2)+1];
                            const V Uk2_re = d[2*(K+N/4)];
                            const V Uk2_im = d[2*(K+N/4)+1];
                            const V Zdk_re = d[2*(K+3*N/4)];
                            const V Zdk_im = d[2*(K+3*N/4)+1];

                            const V w1Zk_re = w1.real()*Zk_re - w1.imag()*Zk_im;
                            const V w1Zk_im = w1.imag()*Zk_re + w1.real()*Zk_im;
                            const V w3Zdk_re = w3.real()*Zdk_re - w3.imag()*Zdk_im;
                            const V w3Zdk_im = w3.imag()*Zdk_re + w3.real()*Zdk_im;

                            d[2*(K)]       = Uk_re  + (w1Zk_re+w3Zdk_re);
                            d[2*(K)+1]     = Uk_im  + (w1Zk_im+w3Zdk_im);

                            d[2*(K+N/2)]   = Uk_re  - (w1Zk_re+w3Zdk_re);
                            d[2*(K+N/2)+1] = Uk_im  - (w1Zk_im+w3Zdk_im);

                            d[2*(K+N/4)]     = Uk2_re + Sign*(-(w1Zk_im-w3Zdk_im));
                            d[2*(K+N/4)+1]   = Uk2_im + Sign*(+(w1Zk_re-w3Zdk_re));

                            d[2*(K+3*N/4)]   = Uk2_re - Sign*(-(w1Zk_im-w3Zdk_im));
                            d[2*(K+3*N/4)+1] = Uk2_im - Sign*(+(w1Zk_re-w3Zdk_re));

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


