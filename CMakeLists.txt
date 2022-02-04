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
        src/PafAlignmentChain.cpp
        src/PafElement.cpp
        )

project(liger2liger)
add_library(liger2liger STATIC ${SOURCES})

# To make sure the library is named Name.so,
# get rid of the "lib" prefix.
set_target_properties(liger2liger PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(liger2liger PROPERTIES DEFINE_SYMBOL "")


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
            )

endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        filter_chimeras_from_paf
        )

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            liger2liger
            Threads::Threads
            )

endforeach()

# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        liger2liger
        DESTINATION ${BINARY_INSTALL_DIR})


