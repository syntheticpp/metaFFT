#pragma once


namespace metaFFT
{
    namespace radix2
    {
        namespace in_place
        {
            // N=2
            template<typename C, template<int, unsigned, class> class butterfly_policy>
            struct radix2<2, C, butterfly_policy>
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

        }
    }
}
