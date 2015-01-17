
#include "bench_helper.h"


#include <fftw3.h>

template<int Sign>
nsec_t fftw_transform(Data* d)
{
    if (d->in == 0) {
        fftwf_complex* in =  (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * d->N);
        d->in = (float*)in;
        d->out = d->in;
        d->plan = (void*)fftwf_plan_dft_1d(d->N, in, in, Sign == -1 ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_PATIENT);
        initData(d);
    }

    const nsec_t t0 = nsec_clock();
    fftwf_execute((fftwf_plan)d->plan);
    const nsec_t dt = nsec_clock() - t0;
    return dt;
}


void fftw_clean(Data* d)
{
    fftwf_free((fftwf_complex*)d->in);
    fftwf_destroy_plan((fftwf_plan)d->plan);
}

