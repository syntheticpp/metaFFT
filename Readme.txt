metaFFT
-------

Template based C++11 Fast-Fourier-Transform implementation.


Idea:
- Completely unroll all loops at compile time with the help of templates.
- Calculate all numerical constants at complile time by using 'constexpr'.
- Use policies for different implementations (complex, Fortran like C, SIMD).


Speed:
Simple  Cooley–Tukey/radix-2 implementation with about 100 lines of code
is 'only' same factors slower than FFTW:

$ ./bin/radix2_sse2_speed
N = 2^ 2 =     4: FFTW/metaFFT =  1.1
N = 2^ 3 =     8: FFTW/metaFFT =  1.2
N = 2^ 4 =    16: FFTW/metaFFT =  1.5
N = 2^ 5 =    32: FFTW/metaFFT =  1.7
N = 2^ 6 =    64: FFTW/metaFFT =  2.0
N = 2^ 7 =   128: FFTW/metaFFT =  2.0
N = 2^ 8 =   256: FFTW/metaFFT =  2.3
N = 2^ 9 =   512: FFTW/metaFFT =  2.8
N = 2^10 =  1024: FFTW/metaFFT =  3.2
N = 2^11 =  2048: FFTW/metaFFT =  3.3
N = 2^12 =  4096: FFTW/metaFFT =  3.4


Build:
- CMake based
- pass -Dlarge=1 to enable large FFTs, this will stress your compiler!


License:
GPL2 with linking exemption.


Links:
http://anthonix.com/ffts
http://nr.com
http://www.drdobbs.com/cpp/a-simple-and-efficient-fft-implementatio/199500857

