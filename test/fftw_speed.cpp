//typedef float float_type;
typedef float float_type;

#include "bench_helper_fftw.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>


int main(int, char*[])
{
    int n = 20;
    int Nmax = 1<<n;

    printf("\nRunning benchmarks until N = 2^%.0f = %i -> arraysize = %ikB\n\n", log2(Nmax), Nmax, 2*Nmax/1024*(int)sizeof(float));
    for (int N = 4; N <= Nmax; N = N<<1) {
        int runs = 1000;
        int nsec = 0;
        float gw = bench(N, "FFTW ", runs, &fftw_transform<-1>, &fftw_clean, nsec);
        if (nsec/1000/1000 > 10) {
            printf("N = 2^%2.0f = %7i: %7i ms, speed FFTW = %4.1f GFLOPS\n", log2(N), N, nsec/1000/1000, gw);
        } else if (nsec/1000 > 10) {
            printf("N = 2^%2.0f = %7i: %7i us, speed FFTW = %4.1f GFLOPS\n", log2(N), N, nsec/1000, gw);
        } else {
            printf("N = 2^%2.0f = %7i: %7i ns, speed FFTW = %4.1f GFLOPS\n", log2(N), N, nsec, gw);
        }
    }

    return 0;
}

