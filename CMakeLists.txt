# CMakeList.txt: PolygonCenterLine 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.17)
project ("PolygonCenterLine")
find_package(CGAL CONFIG REQUIRED COMPONENTS Qt5)
add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)

# 将源代码添加到此项目的可执行文件

file(GLOB header_files "include/*.h")
#message("header_files:" ${header_files})
set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/install" CACHE STRING "" FORCE)
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")

set(CMAKE_CXX_STANDARD 11)
add_executable (PolygonCenterLine "src/PolygonCenterLine.cpp" ${header_files})
target_include_directories(PolygonCenterLine PRIVATE include)
target_link_libraries(PolygonCenterLine PUBLIC CGAL::CGAL PUBLIC CGAL::CGAL_Qt5 PUBLIC)
