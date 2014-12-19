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
                    V* d = (V*)&in[0];
                    V* out = (V*)&o[0];
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
                    V* d = (V*)&in[0];
                    V* out = (V*)&o[0];
                    out[0] = d[0] + d[2*Stride];
                    out[1] = d[1] + d[2*Stride+1];
                    out[2] = d[0] - d[2*Stride];
                    out[3] = d[1] - d[2*Stride+1];
                }
            };

        }
    }
}
