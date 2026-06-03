cmake_minimum_required(VERSION 3.0.0)
project(RESMSearch)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Comprehensive optimization flags for best performance
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -ffast-math -flto -fno-math-errno -finline-functions -foptimize-sibling-calls")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

find_package(PkgConfig REQUIRED)
pkg_check_modules(RNALIB REQUIRED RNAlib2)
find_package(OpenMP)

set(RESM_CORE_SOURCES
    src/buildmap.cpp
    src/align.cpp
    src/msa.cpp
    src/load_data.cpp
)

add_library(resm_core STATIC ${RESM_CORE_SOURCES})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${RNALIB_INCLUDE_DIRS})
link_directories(${RNALIB_LIBRARY_DIRS})

target_include_directories(resm_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${RNALIB_INCLUDE_DIRS}
)
target_compile_options(resm_core PUBLIC ${RNALIB_CFLAGS_OTHER})

target_link_libraries(resm_core PUBLIC ${RNALIB_LDFLAGS} pthread)
if(OpenMP_CXX_FOUND)
    target_link_libraries(resm_core PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(RESM_Search apps/main.cpp)
add_executable(cmapbuild apps/cmapbuild.cpp)
add_executable(cmapsearch apps/cmapsearch.cpp)
add_executable(RESM_Search_align apps/RESM_Search_align.cpp)
add_executable(RESM_Search_alignNCBI apps/RESM_Search_alignNCBI.cpp)
add_executable(BPmap apps/BPmap.cpp)

target_link_libraries(RESM_Search PRIVATE resm_core)
target_link_libraries(cmapbuild PRIVATE resm_core)
target_link_libraries(cmapsearch PRIVATE resm_core)
target_link_libraries(RESM_Search_align PRIVATE resm_core)
target_link_libraries(RESM_Search_alignNCBI PRIVATE resm_core)
target_link_libraries(BPmap PRIVATE resm_core)
