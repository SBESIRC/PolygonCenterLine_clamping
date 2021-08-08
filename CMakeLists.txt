# CMakeList.txt: PolygonCenterLine 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.17)
project ("PolygonCenterLine")
find_package(CGAL CONFIG REQUIRED COMPONENTS Qt5)
find_package(jsoncpp REQUIRED)
add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)

# 将源代码添加到此项目的可执行文件

file(GLOB header_impl "include/impl/*.h")
file(GLOB test_headers "test/*.h")
file(GLOB lib_src "src/*.cpp")
#message("header_files:" ${header_files})
set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/install" CACHE STRING "" FORCE)
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")

set(CMAKE_CXX_STANDARD 11)

add_library(PolygonCenterLine STATIC ${lib_src} ${header_files} ${header_impl})
target_include_directories(PolygonCenterLine PUBLIC include PUBLIC "./libraries")
target_link_libraries(PolygonCenterLine PUBLIC CGAL::CGAL PUBLIC CGAL::CGAL_Qt5 PUBLIC jsoncpp_lib)

# only for testing
add_executable (PolygonCenterLineTest "test/PolygonCenterLineTest.cpp"  ${test_headers})
target_include_directories(PolygonCenterLineTest PRIVATE test)
#target_link_libraries(PolygonCenterLineTest PUBLIC CGAL::CGAL PUBLIC CGAL::CGAL_Qt5)
target_link_libraries(PolygonCenterLineTest PUBLIC PolygonCenterLine)
