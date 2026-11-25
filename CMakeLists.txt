cmake_minimum_required(VERSION 3.16)
project(option_pricing CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_library(pricing src/pricing.cpp)
target_include_directories(pricing PUBLIC include)
target_compile_options(pricing PRIVATE -Wall -Wextra)

add_executable(option_pricing src/main.cpp)
target_link_libraries(option_pricing PRIVATE pricing)

# GoogleTest (downloaded automatically at configure time)
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_executable(pricing_tests tests/test_pricing.cpp)
target_link_libraries(pricing_tests PRIVATE pricing GTest::gtest_main)