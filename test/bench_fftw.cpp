typedef double float_type;

#include "bench_helper_fftw.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>


int main(int, char*[])
{
    int Nmax = 1 << 12;

    printf("\nRunning benchmarks until N = 2^%.0f = %i = %ikB\n\n", log2(Nmax), Nmax, 2*Nmax/1024*(int)sizeof(double));
    for (int N = 4; N <= Nmax; N = N<<1) {
        int runs = 1000;
        int nsec = 0;
        double gw = bench(N, "FFTW ", runs, &fftw_transform<-1>, &fftw_clean, nsec);
        printf("N = 2^%2.0f = %7i: %7i nsec, speed CTFFT = %4.1f GFLOPS\n", log2(N), N, nsec, gw);
    }

    return 0;
}

