project(ExtendedKF)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

set(sources
    src/main.cpp
    src/tools.cpp
    src/tools.h
    src/Kalman.hpp
    src/Kalman.cpp
    src/Tracking2dFusion.hpp
    src/Tracking2dFusion.cpp)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1*/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

add_executable(ExtendedKF ${sources})

find_library(WS_LIB uWS HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib)

target_link_libraries(ExtendedKF z ssl uv "${WS_LIB}")
target_include_directories(ExtendedKF PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
