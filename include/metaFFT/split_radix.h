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
            template<unsigned N, typename C, int Stride, template<int, unsigned, class> class butterfly_policy>
            struct split_radix
            {
                template<int Sign>
                static void calc(C* in, C* out)
                {
                    typedef split_radix<N/2, C, Stride*2, butterfly_policy> even;
                    typedef split_radix<N/4, C, Stride*4, butterfly_policy> odd;

                    even::template calc<Sign>(in, out);
                    odd::template calc<Sign>(in+Stride, out+N/2);
                    odd::template calc<Sign>(in+3*Stride, out+3*N/4);

                    butterfly_policy<Sign, N, C>::loop(out);
                }
            };


            // fft split-radix
            template<unsigned N, typename C, template<int, unsigned, class> class butterfly_policy>
            struct fft
            {
                static void backward(C* in, C* out)
                {
                    split_radix<N, C, 1, butterfly_policy>::template calc<+1>(in, out);
                }

                static void forward(C* in, C* out)
                {
                    split_radix<N, C, 1, butterfly_policy>::template calc<-1>(in, out);
                }
            };



        }
    }
}
