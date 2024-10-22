cmake_minimum_required(VERSION 3.24.0)
project(astrometry-net-lite VERSION 0.95.0 LANGUAGES C CXX)


# Options
option(ASTROMETRY_NET_LITE_BUILD_EXAMPLE "Build the example (default=ON)" ON)


# Output directories
if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
endif()

if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
endif()


# Download dependencies
include(cmake/dependencies.cmake)


# Compile the library
add_subdirectory(src)


# Compile the example
if (ASTROMETRY_NET_LITE_BUILD_EXAMPLE)
    add_subdirectory(example)
endif()
