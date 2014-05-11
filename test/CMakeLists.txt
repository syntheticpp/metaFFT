

find_path(FFTW_INCLUDES fftw3.h PATHS ${fftwdir})
find_library(FFTW_LIBRARIES NAMES fftw3 libfftw3-3 PATHS ${fftwdir})
if(FFTW_INCLUDES AND FFTW_LIBRARIES)
    set(fftw 1)
    message(STATUS "FFTW found in '${fftwdir}'")
    add_definitions(-DHAVE_FFTW)
    include_directories(${FFTW_INCLUDES})
    if(WIN32)
        file(GLOB dlls ${fftwdir}/*.dll)
        foreach(dll ${dlls})
            if(CMAKE_CONFIGURATION_TYPES)
                foreach(type ${CMAKE_CONFIGURATION_TYPES})
                    file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/${type})
                    configure_file(${dll} ${EXECUTABLE_OUTPUT_PATH}/${type} COPYONLY)
                endforeach()
            else()
                configure_file(${dll} ${EXECUTABLE_OUTPUT_PATH}/${type} COPYONLY)
            endif()
        endforeach()
    endif()
else()
    message("FFTW not found: ${FFTW_INCLUDES}, ${FFTW_LIBRARIES}, set fftwdir")
endif()


# huge template processing
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth=5000")
endif()

macro(test name type)
    add_executable(${name} ${name}${type} ${ARGN})
    if(FFTW_LIBRARIES)
        target_link_libraries(${name} ${FFTW_LIBRARIES})
    endif()
    if(UNIX)
        target_link_libraries(${name} -lpthread -ldl)
    endif()
endmacro()

include_directories(${CMAKE_SOURCE_DIR}/include)

file(GLOB headers ${CMAKE_SOURCE_DIR}/include/metaFFT/*.h)
file(GLOB helpers ${CMAKE_SOURCE_DIR}/test/*.h)

list(APPEND headers ${helpers})

test(test_precision .cpp ${headers})

if(fftw)
    test(bench_fftw .cpp ${headers})
    test(compare_fftw .cpp ${headers})
endif()
