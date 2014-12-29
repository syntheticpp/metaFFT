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
        namespace out_place
        {
            // terminate recursion

            // N=1
            template<typename C, int Stride, template<int, unsigned, class, bool> class butterfly_policy>
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
            template<typename C,  int Stride, template<int, unsigned, class, bool> class butterfly_policy>
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




            // N=4 siwtches to split format
            template<typename C,  int Stride, template<int, unsigned, class, bool> class butterfly_policy>
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

                    // data always in interleaved format
                    Complex* o = (Complex*)&out[0];

                    const Complex hp  = { o[2].r + o[3].r, o[2].i + o[3].i };
                    //const Complex hm  = { o[2].r - o[3].r, o[2].i - o[3].i };
                    const Complex ih  = { -o[2].i + o[3].i, o[2].r - o[3].r }; // { -hm.r, hm.r };

                    const Complex t0 = { o[0].r + hp.r, o[0].i + hp.i };
                    const Complex t1 = { o[0].r - hp.r, o[0].i - hp.i };
                    const Complex t2 = { o[1].r - ih.r, o[1].i - ih.i };
                    const Complex t3 = { o[1].r + ih.r, o[1].i + ih.i };

                    if (Stride == 1) {
                        o[0] = t0;
                        o[1] = t1;
                        o[2] = t2;
                        o[3] = t3;
                    } else {
                        // switch to split format
                        metaFFT::avx::Splitted* i = (metaFFT::avx::Splitted*)&out[0];
                        i[0] = { { t0.r, t2.r,  t1.r, t3.r }, { t0.i, t2.i,  t1.i, t3.i } };
                    }
                }
            };


            // N=4 siwtches to split format
            template<typename C,  int Stride, template<int, unsigned, class, bool> class butterfly_policy>
            struct split_radix<8, C, Stride, butterfly_policy>
            {
                typedef typename C::value_type V;

                template<int Sign>
                static void calc(C* in, C* out)
                {
                    constexpr int N = 8;
                    typedef split_radix<N/2, C, Stride*2, butterfly_policy> even;
                    typedef split_radix<N/4, C, Stride*4, butterfly_policy> odd;

                    even::template calc<Sign>(in, out);
                    odd::template calc<Sign>(in+Stride, out+N/2);
                    odd::template calc<Sign>(in+3*Stride, out+3*N/4);

                    constexpr C I(0,1);
                    C o[8];
                    C*const cplx = (C*)&out[0];
                    metaFFT::avx::Splitted*const spl = (metaFFT::avx::Splitted*)&out[0];
                    {

                        C Uk = { spl->r[0], spl->i[0] }; //out[0].real()) + out[2])*I;}
                        C Zk = cplx[4];
                        C Uk2 = { spl->r[2], spl->i[2] }; //out[1]) + out[3])*I;
                        C Zdk = cplx[6];

                        o[0] = Uk  + (Zk+Zdk);
                        o[4] = Uk  - (Zk+Zdk);
                        o[2] = Uk2 - I*(Zk-Zdk);
                        o[6] = Uk2 + I*(Zk-Zdk);
                    }

                    {
                        C Uk = { spl->r[1], spl->i[1] };
                        C Zk = cplx[5];
                        C Uk2 = { spl->r[3], spl->i[3] }; //out[1]) + out[3])*I;
                        C Zdk = cplx[7];

                        const int K = 0;
                        constexpr C w1 = metaFFT::polar<2*(K+0), N, C, Sign == -1>();
                        constexpr C w3 = metaFFT::polar<2*3*(K+3), N, C, Sign == -1>();

                        o[0] = Uk  + (w1*Zk+w3*Zdk);
                        o[4] = Uk  - (w1*Zk+w3*Zdk);
                        o[2] = Uk2 - I*(w1*Zk-w3*Zdk);
                        o[6] = Uk2 + I*(w1*Zk-w3*Zdk);
                    }

                    if (Stride == 1) {
                        for(int i = 0; i< 8; i++)
                            out[i] = o[i];
                    } else {
                        out[0] = o[0].real() + o[1].real()*I;
                        out[1] = o[2].real() + o[3].real()*I;
                        out[2] = o[0].imag() + o[1].imag()*I;
                        out[3] = o[2].imag() + o[3].imag()*I;
                        out[4] = o[4].real() + o[5].real()*I;
                        out[5] = o[6].real() + o[7].real()*I;
                        out[6] = o[4].imag() + o[5].imag()*I;
                        out[7] = o[6].imag() + o[7].imag()*I;
                    }
                }
            };


        }
    }
}
