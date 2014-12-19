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
    namespace radix2
    {
        namespace std_complex
        {

            template<unsigned N, class C>
            struct bit_reverse_policy
            {
                static void bit_reverse(C* data)
                {
                    unsigned j = 1;
                    for (unsigned i = 1; i < N; i++) {
                        if (j > i) {
                            std::swap(data[j - 1], data[i - 1]);
                        }
                        unsigned m = N / 2;
                        while (m >= 2 && j > m) {
                            j -= m;
                            m >>= 1;
                        }
                        j += m;
                    }
                }
            };


            namespace loop
            {
                template<int Sign, unsigned N, class C>
                struct butterfly_policy
                {
                    static void loop(C* d)
                    {
                        constexpr C h = metaFFT::polar<2, N, C, Sign == -1>();
                        C p(1, 0);
                        for (unsigned k = 0; k < N/2; k++, p *= h) {
                            const C o = p * d[k + N/2]; // p * odd element
                            d[k + N/2] = d[k] - o;
                            d[k] += o;
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
                        static void steps(C* d)
                        {
                            const C o = d[K + N/2] * metaFFT::polar<2*K, N, C, Sign == -1>();
                            d[K + N/2] = d[K] - o;
                            d[K] += o;
                            remaining<K+1, End>::steps(d);
                        }

                    };

                    template<unsigned End>
                    struct remaining<End, End> { static void steps(C*) {} };


                    static void loop(C* data)
                    {
                        remaining<0, N/2>::steps(data);
                    }

                };
            }
        }
    }
}


