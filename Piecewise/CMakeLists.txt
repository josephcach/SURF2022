add_library(PWTrilinear ${CMAKE_CURRENT_SOURCE_DIR}/Piecewise_trilinear.cpp)
set_property(TARGET PWTrilinear PROPERTY CXX_STANDARD 20)

set(GTEST_PATH ${CMAKE_SOURCE_DIR}/../googletest)
add_subdirectory(${GTEST_PATH} ${CMAKE_BINARY_DIR}/gtest)

#here you have the tests, now set up for chebyshev
add_executable(Main_PWTrilinear ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
set_property(TARGET Main_PWTrilinear PROPERTY CXX_STANDARD 20)
target_link_libraries(Main_PWTrilinear PWTrilinear)

add_executable(Test_PWTrilinear ${CMAKE_SOURCE_DIR}/Testframe.cpp ${CMAKE_SOURCE_DIR}/Piecewise_trilinear.cpp)
set_property(TARGET Test_PWTrilinear PROPERTY CXX_STANDARD 20)
target_include_directories(Test_PWTrilinear SYSTEM PUBLIC ${GTEST_PATH}/googletest/include)
target_link_libraries(Test_PWTrilinear gtest_main)