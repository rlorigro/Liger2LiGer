######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

set(CMAKE_FIND_DEBUG_MODE TRUE)

cmake_minimum_required(VERSION 3.10)
project(liger2liger VERSION 0.0.0)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")

# Compiler options.
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS )
#set(ASAN_OPTIONS=check_initialization_order=1)
#set(ASAN_OPTIONS=detect_leaks=1)

set(CMAKE_CXX_FLAGS "-fsanitize=address -fexceptions")  # leak checking and address sanitization
#set(CMAKE_CXX_FLAGS "-fexceptions")

add_definitions(-ggdb3 -O0 -Wall)       # Debugging + No optimization
#add_definitions(-O3 -Wall)              # Much optimization

# Definitions needed to eliminate runtime dependency
# on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")


#########################################
# ------------------------------------- #
# -------- SOURCES AND HEADERS -------- #
# ------------------------------------- #
#########################################

# Include header files
include_directories(
        "inc"
)

# Define our shared library sources. NOT test/executables.
set(SOURCES
        src/AlignmentChain.cpp
        src/PafElement.cpp
        src/Bam.cpp
        src/Sam.cpp
        )

project(liger2liger)
add_library(liger2liger STATIC ${SOURCES})

# To make sure the library is named Name.so,
# get rid of the "lib" prefix.
set_target_properties(liger2liger PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(liger2liger PROPERTIES DEFINE_SYMBOL "")

##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ----------------- htslib ----------------- #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

# --------------------

message(STATUS "Fetching htslib")

# https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)

include(FetchContent)

FetchContent_Declare(
        project_htslib
        URL      https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2
)

FetchContent_MakeAvailable(project_htslib)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/htslib/
#        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/os/
#        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/cram/
#        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/htscodecs/
#        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/htscodecs/htscodecs/
        ${CMAKE_SOURCE_DIR}/external/
)

message(STATUS "htslib_SOURCE_DIR: ${project_htslib_SOURCE_DIR}")


set(HTS_DIR _deps/project_htslib-src)

# Do not rebuild if -Ddev=1 is specified in cmake args
if (NOT DEFINED dev OR NOT dev)
    add_custom_target(
            BUILD_HTS
            ALL
            WORKING_DIRECTORY ${HTS_DIR}
            COMMAND pwd
            COMMAND autoreconf -i
            COMMAND ./configure --disable-libcurl --disable-s3 --disable-gcs --without-libdeflate
            COMMAND $(MAKE) print-config
            COMMAND $(MAKE) prefix=${CMAKE_SOURCE_DIR}/external/htslib/ install
    )

endif()

add_library(htslib STATIC IMPORTED)
set_property(TARGET htslib
        PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/external/htslib/lib/libhts.a)

add_dependencies(htslib BUILD_HTS)
add_dependencies(liger2liger htslib)

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")


############################################
# ---------------------------------------- #
# -------- Generating executables -------- #
# ---------------------------------------- #
############################################

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

# -------- TESTS --------

set(TESTS
        )

foreach(FILENAME_PREFIX ${TESTS})
    add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            liger2liger
            Threads::Threads
            htslib
            LibLZMA::LibLZMA
            BZip2::BZip2
            ZLIB::ZLIB
            )

endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        filter_chimeras_from_alignment
        )

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            liger2liger
            Threads::Threads
            htslib
            LibLZMA::LibLZMA
            BZip2::BZip2
            ZLIB::ZLIB

            )

endforeach()

# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        liger2liger
        DESTINATION ${BINARY_INSTALL_DIR})


