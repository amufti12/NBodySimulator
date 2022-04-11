project(NBodySim)

set(NBS_SRCS
    "NBodySim.cpp"
    ${NBS_ROOT_DIR}/thirdparty/imgui/imgui.cpp
    ${NBS_ROOT_DIR}/thirdparty/imgui/imgui_demo.cpp
    ${NBS_ROOT_DIR}/thirdparty/imgui/imgui_draw.cpp
    ${NBS_ROOT_DIR}/thirdparty/imgui/imgui_tables.cpp
    ${NBS_ROOT_DIR}/thirdparty/imgui/imgui_widgets.cpp
    ${NBS_ROOT_DIR}/thirdparty/imgui/backends/imgui_impl_glfw.cpp
    ${NBS_ROOT_DIR}/thirdparty/imgui/backends/imgui_impl_opengl3.cpp
    ${NBS_ROOT_DIR}/thirdparty/implot/implot.cpp
    ${NBS_ROOT_DIR}/thirdparty/implot/implot_items.cpp
    ${NBS_ROOT_DIR}/thirdparty/glad/src/glad.c
)

add_executable(NBodySim ${NBS_SRCS})
target_compile_features(NBodySim PUBLIC cxx_std_20)
target_compile_options(NBodySim PRIVATE ${NBS_CXX_FLAGS})
target_compile_definitions(NBodySim PRIVATE -D_CRT_SECURE_NO_WARNINGS)
target_include_directories(NBodySim PRIVATE 
  ${NBS_ROOT_BUILD_DIR}/config
  ${NBS_ROOT_DIR}/thirdparty/eigen-3.4.0
  ${NBS_ROOT_DIR}/thirdparty/glfw/include
  ${NBS_ROOT_DIR}/thirdparty/imgui
  ${NBS_ROOT_DIR}/thirdparty/implot
  ${NBS_ROOT_DIR}/thirdparty/glad/include
  ${OPENGL_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(NBodySim PRIVATE
  glfw
  ${OpenGL_LIBRARIES}

)
