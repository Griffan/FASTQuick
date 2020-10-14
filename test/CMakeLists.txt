cmake_minimum_required(VERSION 2.8.4)

project(Test)

set(CMAKE_CXX_STANDARD 11)

find_package(GTest REQUIRED)

if (NOT GTEST_FOUND)
    message(FATAL_ERROR "Google Testing framework not found")
endif()

enable_testing()

add_executable(runUnitTests alltests.cpp)

target_link_libraries(runUnitTests GTest::GTest GTest::Main RegionList.cpp)
#add_test(runUnitTests runUnitTests) #included in all tutorials but I don't know what it actually does.
