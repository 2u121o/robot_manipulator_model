cmake_minimum_required(VERSION 3.5)

project(robot_maipulator_model)
set(library_name robot_model)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(include)

set(SOURCES
    src/dynamical_model.cpp
    src/kinematic_model.cpp
    src/inverse_kinematic.cpp
    src/main.cpp
)

set(HEADERS
    include/dynamical_model.h
    include/kinematic_model.h
    include/inverse_kinematic.h
)

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )


add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
add_library(${library_name} SHARED ${SOURCES} ${HEADERS})

#target_link_libraries(${library_name})

install(TARGETS ${library_name}
    LIBRARY DESTINATION /usr/local/lib
    ARCHIVE DESTINATION /usr/local/lib
)

install(DIRECTORY include
    DESTINATION /usr/local
    FILES_MATCHING PATTERN "*.h"
)