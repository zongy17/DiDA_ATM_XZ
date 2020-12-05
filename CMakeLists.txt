cmake_minimum_required(VERSION 3.17)
project(XZ_v11)

set(CMAKE_CXX_STANDARD 11)

#############################################
# 把Armadillo的头文件的目录包括进来（固定的）
INCLUDE_DIRECTORIES("D:/Armadillo/include")
# 把Armadillo lib 目录路径 链接进来
link_directories("D:/Armadillo/lib")
# 链接 Armadillo lib，把那两个lib路径写进来
link_libraries(
        "D:/Armadillo/lib/lapack_win64_MT.lib"
        "D:/Armadillo/lib/blas_win64_MT.lib"
)
#############################################

#add_definitions(-DDEBUG)
add_definitions(-DOUTPUT)
add_definitions(-DPERIODIC)
#add_definitions(-DUPSTREAM)

add_definitions(-Dmod_RK3)

#add_definitions(-DCMAKE_BUILD_TYPE=Release)
#add_definitions(-DCMAKE_BUILD_TYPE=Debug)

add_executable(XZ_v11 main.cpp params.h State.cpp State.h damping.cpp diagnostic.cpp integrate.h integrate.cpp interp.cpp print_func.cpp tendency.cpp)