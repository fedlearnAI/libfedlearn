cmake_minimum_required(VERSION 3.10)
project(testHME VERSION 0.1)

set(CMAKE_CXX_STANDARD 11)

find_package(JNI REQUIRED)
# 加入jni支持
include_directories(${JAVA_INCLUDE_PATH})
include_directories(${JAVA_INCLUDE_PATH2})

include_directories("/usr/include")
include_directories("./include")
include_directories("/usr/include/x86_64-linux-gnu/")

FILE(GLOB SC_FILES "./src/*")

FILE(GLOB NTL_STATIC_LIB "./third_party_lib/libntl/lib/*.a")
message("NTL_STATIC_LIB = ${NTL_STATIC_LIB}")

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
target_link_libraries(Distpaillier ${NTL_STATIC_LIB} m gmp pthread -static-libgcc -static-libstdc++)
set_target_properties(Distpaillier PROPERTIES VERSION ${PROJECT_VERSION})


# build test files
set(i 0)
find_package(GTest REQUIRED)
file(GLOB TEST_FILES "tests/*.cpp")
FOREACH (file ${TEST_FILES})
    add_executable("validate_${i}" ${SC_FILES} ${file})
    target_link_libraries("validate_${i}" ${NTL_STATIC_LIB} ${GTEST_BOTH_LIBRARIES} m gmp pthread gtest -static-libgcc -static-libstdc++)
    math(EXPR i "${i} + 1")
ENDFOREACH ()