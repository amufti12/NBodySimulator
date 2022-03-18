project(NBodySim)

set(NBS_SRCS
    "NBodySim.cpp"
)

add_executable(NBodySim ${NBS_SRCS})
target_compile_features(NBodySim PUBLIC cxx_std_20)
target_compile_options(NBodySim PRIVATE ${NBS_CXX_FLAGS})
target_compile_definitions(NBodySim PRIVATE -D_CRT_SECURE_NO_WARNINGS)
target_include_directories(NBodySim PRIVATE 
  ${NBS_ROOT_BUILD_DIR}/config
  ${CMAKE_CURRENT_SOURCE_DIR}
)
