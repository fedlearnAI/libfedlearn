cmake_minimum_required(VERSION 3.10)
project(libfedlearn VERSION 0.1)

set(CMAKE_CXX_STANDARD 11)

find_package(JNI REQUIRED)
# 加入jni支持
include_directories(${JAVA_INCLUDE_PATH})
include_directories(${JAVA_INCLUDE_PATH2})
include_directories("/usr/include")
include_directories("/usr/include/x86_64-linux-gnu/")
include_directories("./thirdparty/ntl-11.4.3/include")
include_directories("./include")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
message("CMAKE_CXX_COMPILER_ID= ${CMAKE_CXX_COMPILER_ID}")
# STATIC_COMPILE_CMD=""
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
	# set(STATIC_COMPILE_CMD "-static")
	set(STATIC_COMPILE_CMD "")
	set(GMP_STATIC_LIB "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/libgmp/lib/libgmp.dylib")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(STATIC_COMPILE_CMD "-static-libgcc -static-libstdc++")
	set(GMP_STATIC_LIB "gmp")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	set(STATIC_COMPILE_CMD "-static-libgcc -static-libstdc++")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	set(STATIC_COMPILE_CMD "-static-libgcc -static-libstdc++")
endif()


include(ExternalProject)
ExternalProject_Add(
  libntl
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build/libntl/
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ntl-11.4.3/src
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ntl-11.4.3/src/configure "PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/build/libntl/"
  BUILD_COMMAND make "CXXFLAGS=-fPIC"
  BUILD_IN_SOURCE 1
  INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build/libntl/lib
)
ExternalProject_Add_StepTargets(libntl install)
ExternalProject_Get_property(libntl INSTALL_DIR)
message("Source dir of myExtProj = ${INSTALL_DIR}")
message("GMP_STATIC_LIB = ${GMP_STATIC_LIB}")

include_directories("./include")

FILE(GLOB SC_FILES "./src/*")

#set(GMP_STATIC_LIB
#	"${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/libgmp/lib/libgmp.dylib")
set(NTL_STATIC_LIB "${INSTALL_DIR}/libntl.a")
message("1NTL_STATIC_LIB = ${NTL_STATIC_LIB}")


find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

# build shared library
add_library(Distpaillier SHARED ${SC_FILES})
add_dependencies(Distpaillier libntl)
message("GMP_STATIC_LIB=${GMP_STATIC_LIB}")
target_link_libraries(Distpaillier ${NTL_STATIC_LIB} m ${GMP_STATIC_LIB} pthread ${STATIC_COMPILE_CMD})
set_target_properties(Distpaillier PROPERTIES VERSION ${PROJECT_VERSION})


# build test files
include(FetchContent)
FetchContent_Declare(
  googletest
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/googletest
  #URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
set(i 0)
file(GLOB TEST_FILES "./test/*.cpp")
FOREACH (file ${TEST_FILES})
	add_executable("validate_${i}" ${SC_FILES} ${file})
    add_dependencies(validate_0 libntl)
	target_link_libraries("validate_${i}" ${NTL_STATIC_LIB}
		${GTEST_BOTH_LIBRARIES} m ${GMP_STATIC_LIB} pthread gtest ${STATIC_COMPILE_CMD})
    math(EXPR i "${i} + 1")
ENDFOREACH ()

include(GoogleTest)
gtest_discover_tests(validate_1)
