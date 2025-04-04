cmake_minimum_required(VERSION 3.20)

# ============================
# SETUP
# ============================
project(Honeycomb2 CXX)
set(Honeycomb2_Version 0.0.1)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_options(
  -Wall -Wextra
  -fPIC
  -O3
  -mtune=native
  -march=native
  -funroll-loops
  -Wno-unused-parameter -Wno-unused-variable # To be removed later
)

# # -flto causes problems with YAML-CPP

# Check the processor architecture
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|i[3-6]86|amd64")
  # x86-specific compiler flags
  add_compile_options(-mfpmath=sse -msse2)
  message(STATUS "CPU architecture is one of: x86_64|i[3-6]86|amd64")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
  # ARM-specific compiler flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native -mtune=native")

  message(STATUS "CPU architecture is one of: arm64|aarch64")
else()
  message(WARNING "Unrecognized CPU architecture")
endif()

if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH 0)
endif()

if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.27.0")
  cmake_policy(SET CMP0148 OLD)
endif()

# ============================
# MAIN
# ============================
# Generate list of source files
file(GLOB source_files
  src/*.cc
  src/kernels/*.cc
)

# Define target library
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
add_library(Honeycomb2 SHARED ${source_files})

target_include_directories(Honeycomb2 PUBLIC ${PROJECT_SOURCE_DIR}/inc)

# Define libraries to be linked
target_link_libraries(Honeycomb2)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found, enabling OpenMP support.")
  target_link_libraries(Honeycomb2 OpenMP::OpenMP_CXX)
  target_compile_definitions(Honeycomb2 PRIVATE USE_OPENMP)
else()
  message(STATUS "OpenMP not found, building without OpenMP support.")
endif()

# Find PkgConfig
find_package(PkgConfig REQUIRED)

pkg_search_module(EIGEN3 eigen3)

if(NOT EIGEN3_FOUND)
  message(FATAL_ERROR "Could not find Eigen3, aborting...")
else()
  message(STATUS "Found Eigen in ${EIGEN3_INCLUDE_DIRS}")

  # a bit of an hack. The EIGEN3_INCLUDE_DIRS contains the path to Eigen, I need the path to eigen3/ (which contains Eigen/)
  target_include_directories(Honeycomb2 PUBLIC "${EIGEN3_INCLUDE_DIRS}/..")
endif()

set(USER_INCLUDE_DIR "" CACHE PATH "Path to search for cereal library (Should be path to parent folder of cereal/)")

find_path(CEREAL_DIR
  NAMES "cereal/cereal.hpp"
  PATHS ${USER_INCLUDE_DIR}
)

if(CEREAL_DIR)
  message(STATUS "Found cereal in: ${CEREAL_DIR}")
  target_include_directories(Honeycomb2 PUBLIC ${CEREAL_DIR})
else()
  message(FATAL_ERROR "Could not find cereal.")
endif()

# Find yaml-cpp
pkg_search_module(YAML_CPP QUIET yaml-cpp)

if(YAML_CPP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${YAML_CPP_CFLAGS}")
  target_link_libraries(Honeycomb2 ${YAML_CPP_LINK_LIBRARIES})
else(YAML_CPP_FOUND)
  message(FATAL_ERROR "yaml-cpp not found!")
endif(YAML_CPP_FOUND)

enable_testing()
add_subdirectory(tests)

# Configure version to file
configure_file("${PROJECT_SOURCE_DIR}/inc/honeycomb2/version.hpp.in" "${PROJECT_SOURCE_DIR}/inc/honeycomb2/version.hpp")

# Configure configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
set(cflags "-I${prefix}/include")
set(ldflags "-L${prefix}/lib -lHoneycomb2")

configure_file("${PROJECT_SOURCE_DIR}/bin/Honeycomb2-config.in" "${PROJECT_BINARY_DIR}/bin/Honeycomb2-config")

# ============================
# INSTALLATION
# ============================

# Installation of code
install(TARGETS Honeycomb2 RUNTIME DESTINATION bin LIBRARY DESTINATION lib)
install(DIRECTORY inc/honeycomb2 DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/bin/Honeycomb2-config DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
