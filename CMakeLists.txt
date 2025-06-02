cmake_minimum_required(VERSION 3.20)

# ============================
# SETUP
# ============================
set(CMAKE_CXX_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS "")

project(Honeycomb2 CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   message(STATUS "No build type specified. Defaulting to Release.")
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

set(Honeycomb2_Version 0.0.1)

# set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_VECTORIZATION "Enable vectorization flags" ON)
if(ENABLE_VECTORIZATION)
   message(STATUS "Vectorization flags will be applied")
else()
   message(STATUS "Vectorization flags will be disabled")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
   set(HONEYCOMB2_CXX_FLAGS
      -Wall -Wextra
      -fPIC
      -O3
      -funroll-loops
   )
   if(ENABLE_VECTORIZATION)
      list(APPEND HONEYCOMB2_CXX_FLAGS  -mtune=native -march=native)
   endif()
   add_compile_options(${HONEYCOMB2_CXX_FLAGS})

elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
   set(HONEYCOMB2_CXX_FLAGS
      -Wall -Wextra -Wpedantic
      -fPIC
      -O2
      -g
      -funroll-loops
   )
   if(ENABLE_VECTORIZATION)
      list(APPEND HONEYCOMB2_CXX_FLAGS  -mtune=native -march=native)
   endif()
   add_compile_options(${HONEYCOMB2_CXX_FLAGS})

   set(SANITIZER_FLAGS "-fsanitize=address,undefined,alignment,bounds,null,vptr")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG ${SANITIZER_FLAGS}")
   set(HONEYCOMB2_CXX_FLAGS "${HONEYCOMB2_CXX_FLAGS} -DDEBUG ${SANITIZER_FLAGS}")
   link_libraries(${SANITIZER_FLAGS})

else()
   set(HONEYCOMB2_CXX_FLAGS
      -Wall -Wextra -Wpedantic
      -fPIC
      -O0
      -g
      -funroll-loops
   )
   if(ENABLE_VECTORIZATION)
      list(APPEND HONEYCOMB2_CXX_FLAGS  -mtune=native -march=native)
   endif()
   set(SANITIZER_FLAGS "-fsanitize=address,undefined,alignment,bounds,null,vptr")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG ${SANITIZER_FLAGS}")
   set(HONEYCOMB2_CXX_FLAGS "${HONEYCOMB2_CXX_FLAGS} -DDEBUG ${SANITIZER_FLAGS}")
   link_libraries(${SANITIZER_FLAGS})
   add_compile_options(${HONEYCOMB2_CXX_FLAGS})
endif()

# # -flto causes problems with YAML-CPP
if(ENABLE_VECTORIZATION)
   if(CMAKE_BUILD_TYPE STREQUAL "Release")
      if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|i[3-6]86|amd64")
         # x86-specific compiler flags
         list(APPEND HONEYCOMB2_CXX_FLAGS -mfpmath=sse -msse2)
         add_compile_options(-mfpmath=sse -msse2)
         message(STATUS "CPU architecture is one of: x86_64|i[3-6]86|amd64")
      elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
         # ARM-specific compiler flags
         list(APPEND HONEYCOMB2_CXX_FLAGS -mcpu=native -mtune=native)
         add_compile_options(-mcpu=native -mtune=native)

         message(STATUS "CPU architecture is one of: arm64|aarch64")
      else()
         message(WARNING "Unrecognized CPU architecture")
      endif()
   endif()
endif()

if(NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif()

if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.27.0")
   cmake_policy(SET CMP0148 OLD)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Build type specific CXX Flags: ${CMAKE_CXX_FLAGS}")

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
target_include_directories(Honeycomb2 PUBLIC ${PROJECT_SOURCE_DIR}/inc/honeycomb2)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/lib/llvm-19/lib/clang/19/include")

# Define libraries to be linked
target_link_libraries(Honeycomb2)

# set(USER_INCLUDE_DIR "" CACHE PATH "Path to search for cereal library (Should be path to parent folder of cereal/)")

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

string(JOIN " " HONEYCOMB2_CXX_FLAGS_STR ${HONEYCOMB2_CXX_FLAGS})
set(allflags "-I${prefix}/include -std=c++20 ${HONEYCOMB2_CXX_FLAGS_STR}")
set(optflags ${HONEYCOMB2_CXX_FLAGS_STR})
# string(JOIN " " HONEYCOMB2_VEC_FLAGS_STR ${HONEYCOMB2_VEC_FLAGS})
# set(vecflags "${HONEYCOMB2_VEC_FLAGS_STR}")

configure_file("${PROJECT_SOURCE_DIR}/bin/Honeycomb2-config.in" "${PROJECT_BINARY_DIR}/bin/Honeycomb2-config")

# ============================
# INSTALLATION
# ============================

# Installation of code
install(TARGETS Honeycomb2 RUNTIME DESTINATION bin LIBRARY DESTINATION lib)
install(DIRECTORY inc/honeycomb2 DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/bin/Honeycomb2-config DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)


find_package(Doxygen QUIET)

if(DOXYGEN_FOUND)
   message(STATUS "Found Doxygen: ${DOXYGEN_EXECUTABLE}")

   set(DOXYGEN_DIR ${CMAKE_SOURCE_DIR}/docs/doxygen)
   set(DOXYFILE ${DOXYGEN_DIR}/Doxyfile)
   set(DOC_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/docs)

   add_custom_target(Honeycomb2-doc
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE}
      WORKING_DIRECTORY ${DOXYGEN_DIR}
      COMMENT "Generating API documentation with Doxygen"
      VERBATIM
   )
else()
   message(WARNING "Doxygen not found in the system! Cannot generate documentation.")
endif()