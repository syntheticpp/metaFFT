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
        namespace in_place
        {
            // N=2
            template<typename C,template<int, unsigned, class> class butterfly_policy>
            struct radix2<2, C, butterfly_policy>
            {
                template<int Sign>
                static void calc(C* data)
                {
                    static_assert(static_cast<typename C::value_type*>((double*)0) == 0, "double needed");
                    double* d = (double*)&data[0];
                    const __m128d a = _mm_loadu_pd(&d[0]);
                    const __m128d b = _mm_loadu_pd(&d[2]);
                    const __m128d ra = _mm_add_pd(a, b);
                    const __m128d rb = _mm_sub_pd(a, b);
                    _mm_storeu_pd(&d[0], ra);
                    _mm_storeu_pd(&d[2], rb);
                }
            };
        }
    }
}
