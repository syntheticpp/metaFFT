/*  Copyright (C) 2011-2014  Povilas Kanapickas <povilas@radix.lt>

    Distributed under the Boost Software License, Version 1.0.
        (See accompanying file LICENSE_1_0.txt or copy at
            http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef LIBSIMDPP_SIMDPP_DETAIL_INSN_F_SQRT_H
#define LIBSIMDPP_SIMDPP_DETAIL_INSN_F_SQRT_H

#ifndef LIBSIMDPP_SIMD_H
    #error "This file must be included through simd.h"
#endif

#include <cmath>
#include <simdpp/types.h>
#include <simdpp/core/f_rsqrt_e.h>
#include <simdpp/core/f_rsqrt_rh.h>
#include <simdpp/detail/null/foreach.h>
#include <simdpp/detail/null/math.h>

namespace simdpp {
#ifndef SIMDPP_DOXYGEN
namespace SIMDPP_ARCH_NAMESPACE {
#endif
namespace detail {
namespace insn {


SIMDPP_INL float32x4 i_sqrt(const float32x4& a)
{
#if SIMDPP_USE_NULL || SIMDPP_USE_NEON_NO_FLT_SP
    return detail::null::foreach<float32x4>(a, [](float a){ return std::sqrt(a); });
#elif SIMDPP_USE_SSE2
    return _mm_sqrt_ps(a);
#elif SIMPDP_USE_NEON64
    return vsqrtq_f32(a);
#elif SIMDPP_USE_NEON_FLT_SP || SIMDPP_USE_ALTIVEC
    float32x4 x;
    x = rsqrt_e(a);
    x = rsqrt_rh(x, a);
    return mul(a, x);
#endif
}

#if SIMDPP_USE_AVX
SIMDPP_INL float32x8 i_sqrt(const float32x8& a)
{
    return _mm256_sqrt_ps(a);
}
#endif

#if SIMDPP_USE_AVX512
SIMDPP_INL float32<16> i_sqrt(const float32<16>& a)
{
    return _mm512_sqrt_ps(a);
}
#endif

template<unsigned N> SIMDPP_INL
float32<N> i_sqrt(const float32<N>& a)
{
    SIMDPP_VEC_ARRAY_IMPL1(float32<N>, i_sqrt, a);
}

// -----------------------------------------------------------------------------

SIMDPP_INL float64x2 i_sqrt(const float64x2& a)
{
#if SIMDPP_USE_NULL || SIMDPP_USE_NEON32 || SIMDPP_USE_ALTIVEC
    return detail::null::foreach<float64x2>(a, [](double a){ return std::sqrt(a); });
#elif SIMDPP_USE_SSE2
    return _mm_sqrt_pd(a);
#elif SIMDPP_USE_NEON64
    return vsqrtq_f64(a);
#endif
}

#if SIMDPP_USE_AVX
SIMDPP_INL float64x4 i_sqrt(const float64x4& a)
{
    return _mm256_sqrt_pd(a);
}
#endif

#if SIMDPP_USE_AVX512
SIMDPP_INL float64<8> i_sqrt(const float64<8>& a)
{
    return _mm512_sqrt_pd(a);
}
#endif

template<unsigned N> SIMDPP_INL
float64<N> i_sqrt(const float64<N>& a)
{
    SIMDPP_VEC_ARRAY_IMPL1(float64<N>, i_sqrt, a);
}


} // namespace insn
} // namespace detail
#ifndef SIMDPP_DOXYGEN
} // namespace SIMDPP_ARCH_NAMESPACE
#endif
} // namespace simdpp

#endif
