cmake_minimum_required(VERSION 2.8)


set(EXECUTABLE_OUTPUT_PATH   ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH      ${CMAKE_BINARY_DIR}/bin)
set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(metaFFT)

option(large "build large FFTs" OFF)

if(large)
    add_definitions(-DLARGE_FFTS)
    message(STATUS "using large FFTs")
else()
    message(STATUS "pass '-Dlarge=1' to enable large FFTs")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
    add_definitions(-D_USE_MATH_DEFINES)
    #set(MSVC_W_DISABLE "/wd4522 /wd4503")
    #set(MSVC_W_DISABLE "${MSVC_W_DISABLE} /wd4244 /wd4267") # int <-> 64 bit casting
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MSVC_W_DISABLE} /we4013")
else()
    if(UNIX)
        set(fpic -fPIC)
    endif()
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${fpic} -std=c99")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${fpic} -std=c++11 -Wall -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -mfpmath=sse")
    if(NOT avx)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops -funroll-all-loops -flto")
    add_definitions(-DM_PI=3.14159265358979323846264338327)
    add_definitions(-DMETAFFT_HAVE_CONSTEXPR)
endif()


set(include_dir "${CMAKE_SOURCE_DIR}/include")

add_subdirectory(test)
