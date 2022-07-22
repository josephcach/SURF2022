add_library(Piecewise ${CMAKE_CURRENT_SOURCE_DIR}/Piecewise_trilinear.cpp)
set_property(TARGET Piecewise PROPERTY CXX_STANDARD 20)


#here you have the tests, now set up for chebyshev
add_executable(Main_PWTrilinear ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
set_property(TARGET Main_PWTrilinear PROPERTY CXX_STANDARD 20)
target_link_libraries(Main_PWTrilinear Piecewise)


add_executable(Test_PWTrilinear ${CMAKE_CURRENT_SOURCE_DIR}/Testframe.cpp)
set_property(TARGET Test_PWTrilinear PROPERTY CXX_STANDARD 20)
target_link_libraries(Test_PWTrilinear gtest_main Piecewise)