cmake_minimum_required(VERSION 3.10)

project(simulate)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
find_library(TENSORFLOW_LIB tensorflow HINT $ENV{HOME}/libtensorflow2/lib)

execute_process (COMMAND bash -c "echo 'exported directory...'; export LIBRARY_PATH=${LIBRARY_PATH}:~/libtensorflow2/lib; export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:~/libtensorflow2/lib")

set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_C_COMPILER "gcc")
# set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -lasan")

add_executable(${PROJECT_NAME} main.cpp libs/blossom5-v2.05.src/PMduals.o libs/blossom5-v2.05.src/PMexpand.o libs/blossom5-v2.05.src/PMinit.o libs/blossom5-v2.05.src/PMinterface.o libs/blossom5-v2.05.src/PMmain.o libs/blossom5-v2.05.src/PMrepair.o libs/blossom5-v2.05.src/PMshrink.o libs/blossom5-v2.05.src/misc.o libs/blossom5-v2.05.src/MinCost/MinCost.o libs/blossom5-v2.05.src/GEOM/GPMinit.o libs/blossom5-v2.05.src/GEOM/GPMinterface.o libs/blossom5-v2.05.src/GEOM/GPMkdtree.o libs/blossom5-v2.05.src/GEOM/GPMmain.o)
target_link_libraries(${PROJECT_NAME} Threads::Threads)
target_include_directories(${PROJECT_NAME} PRIVATE libs/cppflow-master/include $ENV{HOME}/libtensorflow2/include) #this is the relative path starting from functions.hpp
target_link_libraries (${PROJECT_NAME} "${TENSORFLOW_LIB}")

message("===================Cmake Testing Output===================")
message("CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message("CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")
message("CMAKE_CXX_COMMAND: ${CMAKE_COMMAND}")
message("CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
message("===================        End         ===================")