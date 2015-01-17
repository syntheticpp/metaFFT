#pragma once

#include "nsec_clock.h"

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cmath>


struct Data
{
    int N;
    float* in;
    float* out;
    void* plan;
};


void initData(Data* d)
{
    for(int i = 0; i < d->N; i++) {
        d->in[2*i]   = 0.0;
        d->in[2*i+1] = 0.0;
        d->out[2*i]   = 0.0;
        d->out[2*i+1] = 0.0;
    }
    d->in[2] = 1.0;
}

template<class T>
struct Allocator
{
    T* alloc(int n);
    void free(T*);
};


template<>
struct Allocator<float>
{
    static float* alloc(int n) { return new float[n]; }
    static void free(float* ptr) { delete [] ptr; }
};


template<typename Func, int Sign, class Alloc>
nsec_t metaFFT_transform(Data* d)
{
    if (d->in == 0) {
        d->in = Alloc::alloc(d->N * 2);
        d->out = d->in;
        initData(d);
    }

    const nsec_t t0 = nsec_clock();
    std::complex<float_type>* data = (std::complex<float_type>*)d->in;
    Sign == -1 ? Func::forward(data) : Func::backward(data);
    const nsec_t dt = nsec_clock() - t0;
    return dt;
}


template<typename Func, int Sign, class Alloc>
nsec_t metaFFT_transformOut(Data* d)
{
    if (d->in == 0) {
        d->in = Alloc::alloc(d->N * 2);
        d->out = Alloc::alloc(d->N * 2);
        initData(d);
    }

    const nsec_t t0 = nsec_clock();
    std::complex<float_type>* in = (std::complex<float_type>*)d->in;
    std::complex<float_type>* out = (std::complex<float_type>*)d->out;
    Sign == -1 ? Func::forward(in, out) : Func::backward(in, out);
    const nsec_t dt = nsec_clock() - t0;
    return dt;
}

template<class Alloc>
void metaFFT_clean(Data* d)
{
    if (d->in != d->out)
        Alloc::free(d->out);
    Alloc::free(d->in);
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


#define PI 3.1415926535897932384626433832795028841971693993751058209

double impulse_error(const Data& d, int sign)
{
    long double delta_sum = 0;
    long double sum = 0;

    for(int i = 0; i < d.N; i ++) {
        long double re;
        long double im;
        const long double phi = 2 * PI * (long double)i / (long double)d.N;
        if(sign == -1) {
            re = cosl(phi);
            im = -sinl(phi);
        } else {
            re = cosl(phi);
            im = sinl(phi);
        }
        sum += re * re + im * im;

        re = re - d.out[2*i];
        im = im - d.out[2*i+1];

        delta_sum += re * re + im * im;
    }

    return sqrtl(delta_sum) / sqrtl(sum);
}
