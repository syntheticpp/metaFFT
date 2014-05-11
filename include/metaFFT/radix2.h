#pragma once

namespace metaFFT
{
    namespace radix2
    {
        namespace in_place
        {
            template<unsigned N, typename C, template<int, unsigned, class> class butterfly_policy>
            struct radix2
            {
                template<int Sign>
                static void calc(C* data)
                {
                    typedef radix2<N/2, C, butterfly_policy> recursion;

                    recursion::template calc<Sign>(data);
                    recursion::template calc<Sign>(data + N/2);

                    butterfly_policy<Sign, N, C>::loop(data);
                }
            };


            // terminate recursion
            template<typename C, template<int, unsigned, class> class butterfly_policy>
            struct radix2<1, C, butterfly_policy>
            {
                template<int Sign>
                static void calc(C*) {}
            };



            // fft bitrevere and radix-2 decimate recursively
            template<unsigned N, typename C,
                template<unsigned, class> class bit_reverse_policy,
                template<int, unsigned, class> class butterfly_policy
            >
            struct fft
            {
                static void bit_reverse(C* data)
                {
                    bit_reverse_policy<N, C>::bit_reverse(data);
                }

                static void backward(C* data)
                {
                    bit_reverse(data);
                    radix2<N, C, butterfly_policy>::template calc<+1>(data);
                }

                static void forward(C* data)
                {
                    bit_reverse(data);
                    radix2<N, C, butterfly_policy>::template calc<-1>(data);
                }
            };



        }
    }
}
