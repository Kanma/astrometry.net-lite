# Require at least C++11
if (NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
else()
    if (CMAKE_CXX_STANDARD LESS 17)
        message(FATAL_ERROR "Loading PNM files requires at least C++17")
    endif()
endif()


include_directories(
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>/include
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>/dependencies
)

add_executable(example main.cpp)
target_link_libraries(example PRIVATE astrometry-net-lite)
