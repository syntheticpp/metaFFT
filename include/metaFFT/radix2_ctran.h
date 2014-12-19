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
        namespace ctran
        {

            template<unsigned N, class C>
            struct bit_reverse_policy
            {
                typedef typename C::value_type V;

                static void bit_reverse(C* data)
                {
                    V* d = (V*)&data[0];
                    unsigned j = 1;
                    for (unsigned i = 1; i < 2*N; i += 2) {
                        if (j>i) {
                            std::swap(d[j-1], d[i-1]);
                            std::swap(d[j],   d[i]);
                        }
                        unsigned m = N;
                        while (m>=2 && j>m) {
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
                    typedef typename C::value_type V;

                    static void loop(C* data)
                    {
                        V* d = (V*)&data[0];

                        constexpr V h = Sign * metaFFT::sin<1, N, V>();
                        constexpr V hr =  -2.0 * h * h;
                        constexpr V hi = Sign * metaFFT::sin<2, N, V>();

                        V ht = h;
                        V pr = 1.0;
                        V pi = 0;
                        for (unsigned k = 0; k < N; k += 2) {
                            const V tr = d[k+N]*pr - d[k+N+1]*pi;
                            const V ti = d[k+N]*pi + d[k+N+1]*pr;
                            d[k+N]   = d[k]   - tr;
                            d[k+N+1] = d[k+1] - ti;
                            d[k]     += tr;
                            d[k+1]   += ti;

                            ht = pr;
                            pr += pr*hr - pi*hi;
                            pi += pi*hr + ht*hi;
                        }
                    }
                };
            }


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
                            constexpr C p = metaFFT::polar<K, N, C, Sign == -1>();
                            const V tr = d[K+N]*p.real() - d[K+N+1]*p.imag();
                            const V ti = d[K+N]*p.imag() + d[K+N+1]*p.real();
                            d[K+N]   = d[K]   - tr;
                            d[K+N+1] = d[K+1] - ti;
                            d[K]     += tr;
                            d[K+1]   += ti;

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


