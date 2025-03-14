cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Locate ROOT using root-config
execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --incdir OUTPUT_VARIABLE ROOT_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
message(STATUS "ROOT CXX Flags: ${ROOT_CXX_FLAGS}")
message(STATUS "ROOT Libraries: ${ROOT_LIBS}")

# Add ROOT include directories
message(STATUS "ROOT Include Directory: ${ROOT_INCLUDE_DIR}")
include_directories(${ROOT_INCLUDE_DIR})

# Add project include directories
include_directories(${CMAKE_SOURCE_DIR}/src)

# Define source files
set(SOURCES
        src/CITIROC_Channel.cxx
        app/main.cxx
        src/TransferFunctions.hxx
        src/ILT.hxx
)

# Define the executable
add_executable(my_project ${SOURCES})

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

# Link ROOT libraries
target_link_libraries(my_project ${ROOT_LIBS})