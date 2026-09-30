cmake_minimum_required(VERSION 3.22)
project(genome-variant-toolkit VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_library(genome_variant_core src/version.cpp)
target_include_directories(genome_variant_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genome_variant_core PUBLIC cxx_std_20)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(genome_variant_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(genome-variant-toolkit src/main.cpp)
target_link_libraries(genome-variant-toolkit PRIVATE genome_variant_core)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
