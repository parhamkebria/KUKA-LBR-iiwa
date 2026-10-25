cmake_minimum_required(VERSION 3.10)
project(KUKA CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(kuka KUKA.cpp)
target_include_directories(kuka PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_kuka test_kuka.cpp)
target_link_libraries(test_kuka kuka)

enable_testing()
add_test(NAME kuka_selftest COMMAND test_kuka)
