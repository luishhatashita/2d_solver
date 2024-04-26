# For now I don't know the minumum version of cmake that I will be requiring.
# Otherwise, I will just need to check the mininum required for the other depen-
# dencies.
cmake_minimum_required(VERSION 3.10)

# Create a project
project(2Dsolver VERSION 1.0)

# Make cmake check the C++ version required.
# Set the variable CMAKE_CXX_STANDARD to 11
# and the variable CMAKE_CXX_STANDARD_REQUIRED to True
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -pg")

# Add libraries' subdirectories.
add_subdirectory(parameters)
add_subdirectory(grid)
add_subdirectory(io)
add_subdirectory(malloc)
add_subdirectory(schemes)
add_subdirectory(norms)

# Add main files 
file(GLOB SOURCE_FILES "${CMAKE_CURRENT_LIST_DIR}/*.cpp" "${CMAKE_CURRENT_LIST_DIR}/*.h")

add_executable(2Dsolver 
    ${SOURCE_FILES} 
)

target_link_libraries(2Dsolver PUBLIC parameters grid schemes)

target_include_directories(2Dsolver PUBLIC ${CMAKE_CURRENT_LIST_DIR})
