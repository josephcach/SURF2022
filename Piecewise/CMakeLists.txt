add_library(PWTrilinear ${CMAKE_CURRENT_SOURCE_DIR}/Piecewise_trilinear.cpp)
set_property(TARGET PWTrilinear PROPERTY CXX_STANDARD 20)


#here you have the tests, now set up for chebyshev
add_executable(Main_PWTrilinear ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
set_property(TARGET Main_PWTrilinear PROPERTY CXX_STANDARD 20)
target_link_libraries(Main_PWTrilinear PWTrilinear)