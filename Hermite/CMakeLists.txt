add_library(Hermite ${CMAKE_CURRENT_SOURCE_DIR}/Interpolator_Hermite.cpp)
set_property(TARGET Hermite PROPERTY CXX_STANDARD 20)


#here you have the tests, now set up for chebyshev
add_executable(Test_Hermite ${CMAKE_CURRENT_SOURCE_DIR}/Testframe.cpp)
set_property(TARGET Test_Hermite PROPERTY CXX_STANDARD 20)
target_link_libraries(Test_Hermite gtest_main Hermite)