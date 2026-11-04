# C++ port of the maia_doublets track finding, as a Marlin processor.
# Build with build.sh.
cmake_minimum_required(VERSION 3.15)
project(MaiaDoublets CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG")

find_package(Marlin REQUIRED)
find_package(LCIO REQUIRED)
find_package(nlohmann_json REQUIRED)

add_library(MaiaDoublets SHARED
  src/Calib.cc
  src/Geometry.cc
  src/MiniDoublets.cc
  src/T2s.cc
  src/T4s.cc
  src/T8s.cc
  src/DoubleLayerTrackingProc.cc
)
target_include_directories(MaiaDoublets PRIVATE include ${Marlin_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS})
target_link_libraries(MaiaDoublets ${Marlin_LIBRARIES} ${LCIO_LIBRARIES} nlohmann_json::nlohmann_json)
target_compile_options(MaiaDoublets PRIVATE -Wall -Wextra)

install(TARGETS MaiaDoublets LIBRARY DESTINATION lib)
