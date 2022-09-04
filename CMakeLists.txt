cmake_minimum_required(VERSION 3.20)
project(rotations_analysis_in_manifold)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(EXE_NAME run.exe)
set(EXECUTABLE_OUTPUT_PATH "./")
set(SOURCES
        ./src/main.cpp
)
include_directories(./include)
if (APPLE)
        include_directories(/opt/homebrew/include)
endif()
include_directories(/usr/include)
include_directories(/usr/local/include)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
add_executable(${EXE_NAME} ${SOURCES})
target_link_libraries(
	${EXE_NAME}
	)
