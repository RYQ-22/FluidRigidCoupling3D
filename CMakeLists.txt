cmake_minimum_required(VERSION 3.10)
project(FluidSim VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)
find_package(OpenCV REQUIRED)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

include_directories(include)
include_directories(external/eigen)
include_directories(external/tinyply)
include_directories(${ZLIB_INCLUDE_DIRS})
include_directories(${OpenCV_INCLUDE_DIRS})

set(SOURCE_FILES
    src/main.cpp
    src/fluidsim.cpp
    src/rigidbody.cpp
    src/geometry.cpp
    src/log.cpp
    external/tinyply/tinyply.cpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# link_directories(/path/to/external/library/lib)
# target_link_libraries(${PROJECT_NAME} externalLibraryName) # Link external library
target_link_libraries(${PROJECT_NAME}
    ${ZLIB_LIBRARIES}
    ${OpenCV_LIBS}
    OpenMP::OpenMP_CXX
    /usr/local/lib/libcnpy.so
)
