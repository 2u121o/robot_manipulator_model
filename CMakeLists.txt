cmake_minimum_required(VERSION 3.5)

project(rnea_robot_model)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(include)

set(SOURCES
    src/dynamical_model.cpp
    src/main.cpp
)

set(HEADERS
    include/dynamical_model.h

)

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )


add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})


target_link_libraries(${PROJECT_NAME})
