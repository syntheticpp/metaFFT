/*
 * Copyright (C) 2014 Peter Kümmel. All rights reserved.
 *
 * This file is part of metaFFT, distributed under the GNU GPL v2 with
 * a Linking Exception. For full terms see the included COPYING file.
 */
#pragma once

#include <emmintrin.h>


namespace metaFFT
{
    namespace radix2
    {

        namespace split
        {
            template<typename C,  int Stride, template<int, unsigned, class, bool> class butterfly_policy>
            struct radix2<2, C, Stride, butterfly_policy>
            {
                typedef typename C::value_type V;

                template<int Sign>
                static void calc(C* data)
                {
                    V* d = (V*)&data[0];
                    V tr = d[2];
                    V ti = d[3];
                    d[2] = d[0] - tr;
                    d[3] = d[1] - ti;
                    d[0] += tr;
                    d[1] += ti;
                }
            };

            template<typename C,  int Stride, template<int, unsigned, class, bool> class butterfly_policy>
            struct radix2<4, C, Stride, butterfly_policy>
            {
                typedef typename C::value_type V;

                template<int Sign>
                static void calc(C*)
                {
                }
            };


        }


    }
}
