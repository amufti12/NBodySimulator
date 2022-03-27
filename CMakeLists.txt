cmake_minimum_required(VERSION 3.12.4)
project(NBodySimulator 
        VERSION 0.1 
        DESCRIPTION "A Program which calculates n-bodies of stuff tbd"
        LANGUAGES CXX)

set(NBS_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(NBS_ROOT_BUILD_DIR ${PROJECT_BINARY_DIR})

set(NBS_ARCHIVE_OUTPUT_DIR "${NBS_ROOT_BUILD_DIR}/lib")
set(NBS_LIBRARY_OUTPUT_DIR "${NBS_ROOT_BUILD_DIR}/lib")
set(NBS_RUNTIME_OUTPUT_DIR "${NBS_ROOT_BUILD_DIR}/bin")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${NBS_ARCHIVE_OUTPUT_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${NBS_LIBRARY_OUTPUT_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${NBS_RUNTIME_OUTPUT_DIR})

set(CMAKE_MODULE_PATH ${NBS_ROOT_DIR}/cmake ${CMAKE_MODULE_PATH})

# Exposed variables
option(NBS_ENABLE_ASAN "Enable AddressSanitizer for the builds" OFF)

# Helper things
include(${NBS_ROOT_DIR}/cmake/PlatformDetection.cmake)
include(${NBS_ROOT_DIR}/cmake/CompilerSetup.cmake)

# Versioning variables
set(NBS_VERSION_MAJOR 0)
set(NBS_VERSION_MINOR 1)
set(NBS_VERSION_PATCH 0)
set(NBS_VERSION_STR "0.1.0")

# Configure Common.h.in for general use
configure_file("${NBS_ROOT_DIR}/cmake/Common.h.in"
               "${NBS_ROOT_BUILD_DIR}/config/Common.h")

add_subdirectory(thirdparty/eigen-3.4.0)

add_subdirectory(src)
