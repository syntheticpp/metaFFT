/*
 * Copyright (C) 2014 Peter Kümmel. All rights reserved.
 *
 * This file is part of metaFFT, distributed under the GNU GPL v2 with
 * a Linking Exception. For full terms see the included COPYING file.
 */
#pragma once

#include <cmath>
#include <complex>

#ifndef METAFFT_HAVE_CONSTEXPR
#define constexpr const
#endif

namespace metaFFT
{


#ifdef METAFFT_HAVE_CONSTEXPR

    template<unsigned A, unsigned B, class T>
    constexpr T sin()
    {
        static_assert(::sin(0) == 0, "constexpr sin() needed");
        return ::sin((T)A*M_PI/B);
    }

    template<unsigned A, unsigned B, class T>
    constexpr T cos()
    {
        static_assert(::cos(0) == 1, "constexpr cos() needed");
        return ::cos((T)A*M_PI/B);
    }

    template<unsigned A, unsigned B, class T, bool Conj = false>
    constexpr T polar()
    {
        static_assert(::sin(0) == 0, "constexpr sin() needed");
        static_assert(::cos(0) == 1, "constexpr cos() needed");
        typedef typename T::value_type V;
        return T(::cos((V)A*M_PI/B), (Conj ? -1 : 1) * ::sin((V)A*M_PI/B));
    }

#else

    template<unsigned A, unsigned B, class T>
    struct Series
    {
        template<unsigned M, unsigned I>
        struct Sin {
            static constexpr T value() { return 1.0 - (A*M_PI/B)*(A*M_PI/B)/M/(M+1) * Sin<M+2, I>::value(); }
        };

        template<unsigned I>
        struct Sin<I, I> {
            static constexpr T value() { return 1.0; }
        };
    };

    template<unsigned A, unsigned B, class T>
    constexpr T sin() { return Series<A, B, T>::template Sin<2, 34>::value() * (A*M_PI/B); }

    template<unsigned A, unsigned B, class T>
    constexpr T cos() { return Series<A, B, T>::template Sin<1, 33>::value(); }

    template<unsigned A, unsigned B, class T, bool Conj = false>
    constexpr T polar()
    {
        typedef typename T::value_type V;
        return T(cos<A, B, V>(), (Conj ? -1 : 1) *  sin<A, B, V>());
    }

#endif

}
