
#include "bench_helper.h"


#include <fftw3.h>

template<int Sign>
nsec_t fftw_transform(Data* d)
{
    if (d->in == 0) {
        fftw_complex* in =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * d->N);
        d->in = (double*)in;
        d->out = d->in;
        d->plan = (void*)fftw_plan_dft_1d(d->N, in, in, Sign == -1 ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_PATIENT);
        initData(d);
    }

    const nsec_t t0 = nsec_clock();
    fftw_execute((fftw_plan)d->plan);
    const nsec_t dt = nsec_clock() - t0;
    return dt;
}


void fftw_clean(Data* d)
{
    fftw_free((fftw_complex*)d->in);
    fftw_destroy_plan((fftw_plan)d->plan);
}

