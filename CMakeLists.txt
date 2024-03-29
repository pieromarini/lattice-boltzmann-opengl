cmake_minimum_required(VERSION 3.1)

project(lattice_fluids)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -g -Wall -O3")

set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# File Gathering
file(GLOB_RECURSE HDRS
	${SRC_DIR}/*.h
	${SRC_DIR}/*.hpp
)
file(GLOB_RECURSE SRCS
	${SRC_DIR}/*.cpp
	${SRC_DIR}/*.c
)

add_executable(${PROJECT_NAME} ${HDRS} ${SRCS})

# Use Modern OpenGL
set(OpenGL_GL_PREFERENCE "GLVND")

include(FindPkgConfig)

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

add_subdirectory(${LIB_DIR}/glfw)

find_package(OpenGL REQUIRED)
if(NOT OPENGL_FOUND)
	message(FATAL_ERROR "OpenGL not found!")
endif(NOT OPENGL_FOUND)

find_package(GLEW REQUIRED)
if(NOT GLEW_FOUND)
	message(FATAL_ERROR "GLEW not found!")
endif(NOT GLEW_FOUND)

# Define the include DIRs
include_directories( ${PROJECT_NAME}
	${lattice_fluids_SOURCE_DIR}/src
	${OPENGL_INCLUDE_DIRS}
	${GLEW_INCLUDE_DIRS}
	${GLFW_INCLUDE_DIRS}
)

# Define the link libraries
target_link_libraries( ${PROJECT_NAME}
	glfw
	${OPENGL_LIBRARIES}
	${GLEW_LIBRARIES}
)

# So we can use resources/ folder from build dir.
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${lattice_fluids_SOURCE_DIR}/resources ${lattice_fluids_BINARY_DIR}/resources RESULT_VARIABLE exitcode)
