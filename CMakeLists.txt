cmake_minimum_required(VERSION 3.20)

project(Fractal_Image_Compression)

set(CMAKE_CXX_FLAGS "-Wall -g")

add_executable(enc src/enc.cpp src/imgio.cpp src/indexing.cpp src/coding.cpp src/miscell.cpp src/split.cpp src/nns.cpp)
add_executable(dec src/dec.cpp src/imgio.cpp src/miscell.cpp)
add_executable(eval src/eval.cpp src/imgio.cpp src/miscell.cpp)

target_compile_features(enc PUBLIC cxx_std_17)
target_compile_features(dec PUBLIC cxx_std_17)
target_compile_features(eval PUBLIC cxx_std_17)
