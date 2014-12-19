typedef double float_type;

#include "bench_helper.h"

#ifdef HAVE_FFTW
#include "bench_helper_fftw.h"
#endif

#include <metaFFT/split_radix.h>
#include <metaFFT/split_radix_ctran.h>
#include <metaFFT/split_radix_complex.h>
#include <metaFFT/split_radix_kernel.h>

template<int Sign, typename Transform, typename Cleanup>
void test(int N, const char* name, Transform transform, Cleanup cleanup)
{
    Data d;
    d.N = N;
    d.in = 0;
    d.out = 0;

    transform(&d);

    printf(" %3d  | %5d | %1.2e | %s\n", Sign, N, impulse_error(d, Sign), name);

    cleanup(&d);
}


template<int Sign, typename Func>
void test(int N, const char* name)
{
     test<Sign>(N, name, metaFFT_transformOut<Func, Sign, Allocator<float_type>>, &metaFFT_clean<Allocator<float_type>>);
}


template<int Sign, int N>
void test()
{
#ifdef HAVE_FFTW
    test<Sign>(N, "FFTW ", &fftw_transform<Sign>, &fftw_clean);
#endif

    using metaFFT::split_radix::out_place::fft;


    {
        using namespace metaFFT::split_radix::std_complex;
        using namespace metaFFT::split_radix::std_complex::loop;
        test<Sign, fft<N, std::complex<float_type>, butterfly_policy>>(N, "std+loop      ");
    }


    {
        using namespace metaFFT::split_radix::std_complex;
        using namespace metaFFT::split_radix::std_complex::unrolled_loop;
        test<Sign, fft<N, std::complex<float_type>, butterfly_policy>>(N, "std+unrolled   ");
    }


    {
        using namespace metaFFT::split_radix::ctran;
        using namespace metaFFT::split_radix::ctran::unrolled_loop;
        test<Sign, fft<N, std::complex<float_type>, butterfly_policy>>(N,"ctran+unrolled  ");
    }

    printf("\n");
}


template<int N>
void test()
{
    test<+1, N>();
    test<-1, N>();
}


int main(int, char *[])
{
    printf(" Sign |      Size |     L2 Error\n");
    printf("------+-----------+-------------\n");

    test<   2>();
    test<   4>();
    test<   8>();
    test<  16>();
    test<  32>();
    test<  64>();
    test< 128>();
    test< 256>();

#if 0 || defined(LARGE_FFTS)
    test< 512>();
    test<1024>();
    test<2048>();
    //test<2048*2>();
#endif

    return 0;
}
