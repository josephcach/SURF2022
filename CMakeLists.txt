cmake_minimum_required(VERSION 3.10)
project(SURF2022 CXX)

set(default_build_type "Release")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++20 -Wall -Wextra -Wno-sign-compare -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -ggdb -O0 -DDEBUG") 


#add gtest. This is required.
set(GTEST_PATH ${CMAKE_SOURCE_DIR}/../googletest)
add_subdirectory(${GTEST_PATH} ${CMAKE_BINARY_DIR}/gtest)

#add the hermite and piecewise stuff
add_subdirectory(${CMAKE_SOURCE_DIR}/Hermite)
add_subdirectory(${CMAKE_SOURCE_DIR}/Piecewise)


#enable_testing()

add_executable(SURF2022_trilinear ${CMAKE_SOURCE_DIR}/main.cpp ${CMAKE_SOURCE_DIR}/Interpolator_trilinear.cpp)
set_property(TARGET SURF2022_trilinear PROPERTY CXX_STANDARD 20)
add_executable(SURF2022_cheb ${CMAKE_SOURCE_DIR}/main_cheb.cpp ${CMAKE_SOURCE_DIR}/Interpolator_chebyshev.cpp)
set_property(TARGET SURF2022_cheb PROPERTY CXX_STANDARD 20)


#here you have the tests, now set up for chebyshev
add_executable(Test_cheb ${CMAKE_SOURCE_DIR}/Testframe_interpolator_chebyshev.cpp ${CMAKE_SOURCE_DIR}/Interpolator_chebyshev.cpp)
set_property(TARGET Test_cheb PROPERTY CXX_STANDARD 20)
target_include_directories(Test_cheb SYSTEM PUBLIC ${GTEST_PATH}/googletest/include)
target_link_libraries(Test_cheb gtest_main)

#add_test(Test_cheb Test_cheb)