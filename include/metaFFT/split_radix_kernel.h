/*
 * Copyright (C) 2014 Peter Kümmel. All rights reserved.
 *
 * This file is part of metaFFT, distributed under the GNU GPL v2 with
 * a Linking Exception. For full terms see the included COPYING file.
 */
#pragma once

namespace metaFFT
{
    namespace split_radix
    {
        namespace out_place
        {
            // terminate recursion

            // N=1
            template<typename C, int Stride, template<int, unsigned, class> class butterfly_policy>
            struct split_radix<1, C, Stride, butterfly_policy>
            {
                typedef typename C::value_type V;

                template<int Sign>
                static void calc(C* in, C* o)
                {
                    V*const d = (V*)&in[0];
                    V*const out = (V*)&o[0];
                    out[0] = d[0];
                    out[1] = d[1];
                }
            };


            // N=2
            template<typename C,  int Stride, template<int, unsigned, class> class butterfly_policy>
            struct split_radix<2, C, Stride, butterfly_policy>
            {
                typedef typename C::value_type V;

                template<int Sign>
                static void calc(C* in, C* o)
                {
                    V*const d = (V*)&in[0];
                    V*const out = (V*)&o[0];
                    out[0] = d[0] + d[2*Stride];
                    out[1] = d[1] + d[2*Stride+1];
                    out[2] = d[0] - d[2*Stride];
                    out[3] = d[1] - d[2*Stride+1];
                }
            };

            struct Complex {
                double r;
                double i;
            };

            /*
            // N=4
            template<typename C,  int Stride, template<int, unsigned, class> class butterfly_policy>
            struct split_radix<4, C, Stride, butterfly_policy>
            {
                typedef typename C::value_type V;

                template<int Sign>
                static void calc(C* in, C* out)
                {
                    constexpr int N = 4;
                    typedef split_radix<N/2, C, Stride*2, butterfly_policy> even;
                    typedef split_radix<N/4, C, Stride*4, butterfly_policy> odd;

                    even::template calc<Sign>(in, out);
                    odd::template calc<Sign>(in+Stride, out+N/2);
                    odd::template calc<Sign>(in+3*Stride, out+3*N/4);

                    Complex* o = (Complex*)&out[0];

                    const Complex hp  = { o[2].r + o[3].r, o[2].i + o[3].i };
                    //const Complex hm  = { o[2].r - o[3].r, o[2].i - o[3].i };
                    const Complex ih  = { -o[2].i + o[3].i, o[2].r - o[3].r }; // { -hm.r, hm.r };

                    const Complex t0 = { o[0].r + hp.r, o[0].i + hp.i };
                    const Complex t1 = { o[0].r - hp.r, o[0].i - hp.i };
                    const Complex t2 = { o[1].r - ih.r, o[1].i - ih.i };
                    const Complex t3 = { o[1].r + ih.r, o[1].i + ih.i };

                    constexpr bool log2stride = log2(Stride) != 0;
                    if (log2stride) {
                        o[0] = { t0.r, t2.r };
                        o[1] = { t1.r, t3.r };
                        o[2] = { t0.i, t2.i };
                        o[3] = { t1.i, t3.i };
                    } else {
                        o[0] = t0;
                        o[1] = t1;
                        o[2] = t2;
                        o[3] = t3;
                    }
                }
            };
            */


        }
    }
}
