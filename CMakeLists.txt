cmake_minimum_required(VERSION 3.5)

project(rnea_robot_model)

include_directories(include)

set(SOURCES
    src/dynamical_model.cpp
    src/link.cpp
    src/main.cpp

)

set(HEADERS
    include/dynamical_model.h
    include/link.h

)

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )


add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})


target_link_libraries(${PROJECT_NAME})
