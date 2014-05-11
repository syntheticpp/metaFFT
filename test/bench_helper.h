#pragma once

#include "nsec_clock.h"

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cmath>


struct Data
{
    int N;
    double* in;
    double* out;
    void* plan;
};


void initData(Data* d)
{
    for(int i = 0; i < d->N; i++) {
        d->in[2*i]   = 0.0;
        d->in[2*i+1] = 0.0;
    }
    d->in[2] = 1.0;
}



template<typename Func, int Sign>
nsec_t metaFFT_transform(Data* d)
{
    if (d->in == 0) {
        d->in = new double[d->N * 2];
        d->out = d->in;
        initData(d);
    }

    const nsec_t t0 = nsec_clock();
    std::complex<float_type>* data = (std::complex<float_type>*)d->in;
    Sign == -1 ? Func::forward(data) : Func::backward(data);
    const nsec_t dt = nsec_clock() - t0;
    return dt;
}


void metaFFT_clean(Data* d)
{
    delete [] d->in;
}


double gigaComplexFLOPS(int N, nsec_t nsec)
{
    if (nsec == 0) {
        return 0;
    }
    return 5.0 * N * log2(N) / nsec;
}


template<typename Transform, typename Cleanup>
double bench(int N, const char* msg, int runs, Transform transform, Cleanup cleanup, int& nsec)
{
    Data d;
    d.N = N;
    d.in = 0;
    d.out = 0;

    nsec_t dt_min = 1<<30;
    int updated = 0;
    for (int i = 0; i < runs; i++) {
        nsec_t dt = transform(&d);
        if (dt_min > dt) {
            dt_min = dt;
            updated++;
        }
    }

    // check if average should be used
    if (updated < runs/2) {
        const nsec_t t0_avg = nsec_clock();
        for (int i = 0; i < runs; i++) {
            transform(&d);
        }
        nsec_t tmin_avg = (nsec_clock() - t0_avg ) / runs;
        dt_min = tmin_avg;
    }
    cleanup(&d);

    double gigaflops = gigaComplexFLOPS(N, dt_min);
    nsec = (int)dt_min;
    (void) msg;
    //printf("N=2^%.0f=%i: %s %.2f Gigaflops, %i usec\n", log2(N), N, msg, gigaflops, usec);
    return gigaflops;
}

