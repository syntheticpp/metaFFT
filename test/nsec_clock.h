#pragma once


#ifdef _WIN32
#include "windows.h"
typedef LONGLONG nsec_t;
#else
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
typedef long nsec_t;
#endif

#include <stdio.h>



inline nsec_t nsec_clock()
{
    nsec_t nsecs = 0;

#ifdef _WIN32
    LARGE_INTEGER large;
    QueryPerformanceCounter(&large);
    nsecs = large.QuadPart;
    nsecs *= 1000*1000*1000;
    QueryPerformanceFrequency(&large);
    nsecs /= large.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    nsecs = 1000*1000*1000 * ts.tv_sec;
    nsecs += ts.tv_nsec;
#endif

    return nsecs;
}
