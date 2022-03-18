include(${NBS_ROOT_DIR}/cmake/PlatformDetection.cmake)

# Comiler detection
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(NBS_CXX_COMPILER "Clang")
  set(NBS_CLANG 1)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  set(NBS_CXX_COMPILER "Clang")
  set(NBS_CLANG 1)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(NBS_CXX_COMILER "GNU")
  set(NBS_GCC 1)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set(NBS_CXX_COMPILER "MSVC")
  set(NBS_MSVC 1)
else()
  message(FATAL_ERROR "Building NBS is only supported with the Clang, AppleClang, GNU, and MSVC comiplers.")
endif()

# General CXX flags
set(NBS_CXX_BASE_FLAGS)
if(${NBS_MSVC})
  # Disabled warnings for MSVC take the form /wd####
  # Reference: https://docs.microsoft.com/en-us/cpp/build/reference/compiler-option-warning-level?view=msvc-170
  set(NBS_CXX_DISABLED_WARNINGS)
  set(NBS_CXX_BASE_FLAGS /W4 ${NBS_CXX_DISABLED_WARNINGS} /WX)
else()
  # Disabled warnings for any variant of clang/gcc take the form -Wno-blah
  # Reference: https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html
  set(NBS_CXX_DISABLED_WARNINGS)
  set(NBS_CXX_BASE_FLAGS -Wall -Wextra -pedantic ${NBS_CXX_DISABLED_WARNINGS} -Werror)
endif()

# Handle ASAN flags if requested
if(${NBS_ENABLE_ASAN})
  if(${NBS_MSVC})
    # If the MSVC version is 16.9 or later
    if("${MSVC_VERSION}" VERSION_GREATER_EQUAL "1929")
      set(NBS_CXX_BASE_FLAGS ${NBS_CXX_BASE_FLAGS} /fsanitize=address)
    else()
      message(WARNING "To use ASAN, please upgrade your MSVC compiler to version 16.9 or later.")
    endif()
  else()
    # Going to assume someone is using a version of GCC/Clang which has ASAN capabilities...
    set(NBS_CXX_BASE_FLAGS ${NBS_CXX_BASE_FLAGS} -fsanitize=address -fno-omit-frame-pointer)
  endif()
endif()

set(NBS_CXX_FLAGS ${NBS_CXX_BASE_FLAGS})
