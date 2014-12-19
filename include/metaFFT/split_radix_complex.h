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
        namespace std_complex
        {

            namespace loop
            {
                template<int Sign, unsigned N, class C>
                struct butterfly_policy
                {
                    static void loop(C* out)
                    {
                        constexpr C i(0, 1.0);
                        constexpr C dir(-Sign, 0);
                        {
                            const C Uk  = out[0];
                            const C Zk  = out[N/2];
                            const C Uk2 = out[N/4];
                            const C Zdk = out[3*N/4];
                            out[0]     = Uk  +   (Zk+Zdk);
                            out[N/2]   = Uk  -   (Zk+Zdk);
                            out[N/4]   = Uk2 - dir*i*(Zk-Zdk);
                            out[3*N/4] = Uk2 + dir*i*(Zk-Zdk);
                        }

                        constexpr C h1 = metaFFT::polar<2*1, N, C, Sign == -1>();
                        constexpr C h3 = metaFFT::polar<2*3, N, C, Sign == -1>();
                        C w1 = C(1, 0) * h1;
                        C w3 = C(1, 0) * h3;
                        for (unsigned k = 1; k < N/4; k++, w1 *= h1, w3 *= h3 ) {
                            const C Uk  = out[k];
                            const C Zk  = out[k+N/2];
                            const C Uk2 = out[k+N/4];
                            const C Zdk = out[k+3*N/4];

                            const C w1Zk = w1*Zk;
                            const C w3Zdk = w3*Zdk;
                            out[k]       = Uk  +       (w1Zk+w3Zdk);
                            out[k+N/2]   = Uk  -       (w1Zk+w3Zdk);
                            out[k+N/4]   = Uk2 - dir*i*(w1Zk-w3Zdk);
                            out[k+3*N/4] = Uk2 + dir*i*(w1Zk-w3Zdk);
                        }
                    };
                };
            }



            namespace unrolled_loop
            {
                template<int Sign, unsigned N, class C>
                struct butterfly_policy
                {
                    template<unsigned K, unsigned End>
                    struct remaining
                    {
                        static void steps(C* out)
                        {
                            constexpr C i(0, 1.0);
                            constexpr C dir(-Sign, 0);

                            constexpr C w1 = metaFFT::polar<2*K, N, C, Sign == -1>();
                            constexpr C w3 = metaFFT::polar<2*3*K, N, C, Sign == -1>();
                            const C Uk  = out[K];
                            const C Zk  = out[K+N/2];
                            const C Uk2 = out[K+N/4];
                            const C Zdk = out[K+3*N/4];

                            const C w1Zk = w1*Zk;
                            const C w3Zdk = w3*Zdk;
                            out[K]       = Uk  +       (w1Zk+w3Zdk);
                            out[K+N/2]   = Uk  -       (w1Zk+w3Zdk);
                            out[K+N/4]   = Uk2 - dir*i*(w1Zk-w3Zdk);
                            out[K+3*N/4] = Uk2 + dir*i*(w1Zk-w3Zdk);
            
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


