cmake_minimum_required(VERSION 3.16)

# Sets project name
project(GoingWithTheFlow LANGUAGES CXX C)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# Sets C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Shared include paths
include_directories(src)
include_directories(libs)

# ==============================================================================
# Main Application with Viewer (requires Qt6 + GLEW)
# ==============================================================================
find_package(Qt6 QUIET COMPONENTS Core Concurrent Xml Widgets OpenGL OpenGLWidgets Gui)

if (Qt6_FOUND)
    add_executable(${PROJECT_NAME}
        # Simulation core
        src/integrator.h src/integrator.cpp
        src/momentum.h src/momentum.cpp
        src/inertia.h src/inertia.cpp
        src/lift_and_drag.cpp
        src/mesh_loader.h
        src/se3.h

        # Viewer
        src/main.cpp
        src/simviewerwidget.h src/simviewerwidget.cpp
        src/graphics/camera.h src/graphics/camera.cpp
        src/graphics/shader.h src/graphics/shader.cpp
        src/graphics/shape.h src/graphics/shape.cpp
        src/graphics/graphicsdebug.h src/graphics/graphicsdebug.cpp

        # Debug widgets (keep for now)
        src/bodydebugwidget.h src/bodydebugwidget.cpp
        src/fluiddebugwidget.h src/fluiddebugwidget.cpp
        src/inertiadebugwidget.h src/inertiadebugwidget.cpp

        # Shader resources
        resources/shaders/shaders.qrc
    )

    # GLEW
    add_library(StaticGLEW STATIC glew/src/glew.c)
    include_directories(${PROJECT_NAME} PRIVATE glew/include)

    target_link_libraries(${PROJECT_NAME} PRIVATE
        Qt::Concurrent
        Qt::Core
        Qt::Gui
        Qt::OpenGL
        Qt::OpenGLWidgets
        Qt::Widgets
        Qt::Xml
        StaticGLEW
    )

    target_include_directories(${PROJECT_NAME} PRIVATE
        Eigen
        util
    )

    # GLEW Windows support
    if (WIN32)
      add_compile_definitions(GLEW_STATIC)
      target_link_libraries(${PROJECT_NAME} PRIVATE
        opengl32
        glu32
      )
    endif()
else()
    message(STATUS "Qt6 not found — only building test targets and run_sim")
endif()

# Silence platform warnings
if (MSVC OR MSYS OR MINGW)
  set(CMAKE_CXX_FLAGS "-Wno-volatile")
endif()
if (APPLE)
  set(CMAKE_CXX_FLAGS "-Wno-deprecated-volatile")
endif()
