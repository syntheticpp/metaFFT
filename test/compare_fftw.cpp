typedef double float_type;

#include "bench_helper.h"
#include "bench_helper_fftw.h"

#include <metaFFT/radix2.h>
#include <metaFFT/radix2_kernel.h>
#include <metaFFT/radix2_complex.h>
#include <metaFFT/radix2_ctran.h>


template<typename Func>
void compare(int N, const char* name, int runs, double fftw_flops)
{
    int nsec;
    double ct = bench(N, name, runs, metaFFT_transform<Func, -1>, &metaFFT_clean, nsec);
    printf("N = 2^%2.0f = %5i: %5i nsec, speed %s = %5.2f GFLOPS  FFTW/metaFFT = %4.1f\n", log2(N), N, nsec, name, ct, fftw_flops / ct);
}


template<int N>
void compare()
{
    const int runs = 10*1000;

    int nsec;
    const double fftw_flops = bench(N, "FFTW ", runs, &fftw_transform<-1>, &fftw_clean, nsec);

    using metaFFT::radix2::in_place::fft;

    {
        using namespace metaFFT::radix2::std_complex;
        using namespace metaFFT::radix2::std_complex::unrolled_loop;
        compare<fft<N, std::complex<float_type>, bit_reverse_policy, butterfly_policy>>(N, "std+unrolled  ", runs, fftw_flops);
    }

    {
        using namespace metaFFT::radix2::std_complex;
        using namespace metaFFT::radix2::std_complex::loop;
        compare<fft<N, std::complex<float_type>, bit_reverse_policy, butterfly_policy>>(N, "std+loop      ", runs, fftw_flops);
    }

    {
        using namespace metaFFT::radix2::ctran;
        using namespace metaFFT::radix2::ctran::unrolled_loop;
        compare<fft<N, std::complex<float_type>, bit_reverse_policy, butterfly_policy>>(N, "ctran+unrolled", runs, fftw_flops);
    }

    {
        using namespace metaFFT::radix2::ctran;
        using namespace metaFFT::radix2::ctran::loop;
        compare<fft<N, std::complex<float_type>, bit_reverse_policy, butterfly_policy>>(N, "ctran+loop    ", runs, fftw_flops);
    }

    printf("\n");
}



int main(int, char *[])
{
    compare<  4>();
    compare<  8>();
    compare< 16>();
    compare< 32>();
    compare< 64>();
    compare<128>();
    compare<256>();

#if 0
    compare<512>();
    compare<1024>();
    compare<2048>();
#endif

    return 0;
}

