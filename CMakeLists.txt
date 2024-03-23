cmake_minimum_required(VERSION 3.10)
project(FluidSim VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(include)
include_directories(external/eigen)

set(SOURCE_FILES
    src/main.cpp
    src/fluidsim.cpp    
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# link_directories(/path/to/external/library/lib)
# target_link_libraries(${PROJECT_NAME} externalLibraryName) # Link external library
