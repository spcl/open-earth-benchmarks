cmake_minimum_required(VERSION 3.1)

project(OpenEarthBenchmarks VERSION 1.0
        DESCRIPTION "Open Earth Benchmarks"
        LANGUAGES CXX)

include_directories(include)
include_directories(kernels)

add_executable(laplace src/laplace.cpp)
add_executable(hdiffsa src/hdiffsa.cpp)
add_executable(hdiffsmag src/hdiffsmag.cpp)
add_executable(hadvuv src/hadvuv.cpp)
add_executable(hadvuv5th src/hadvuv5th.cpp)
add_executable(fastwavesuv src/fastwavesuv.cpp)
add_executable(uvbke src/uvbke.cpp)
add_executable(p_grad_c src/p_grad_c.cpp)
add_executable(nh_p_grad src/nh_p_grad.cpp)
add_executable(fvtp2d src/fvtp2d.cpp)
