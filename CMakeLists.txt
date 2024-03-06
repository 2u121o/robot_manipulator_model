cmake_minimum_required(VERSION 3.5)

project(robot_maipulator_model)
set(library_name robot_model)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_CXX_FLAGS_DEBUG "-g")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3")

option(ENABLE_TESTING "Enable google test" OFF)
if(ENABLE_TESTING MATCHES ON)
    add_subdirectory(test)
endif()

include_directories(include)
include_directories(/usr/include/jsoncpp/)

set(SOURCES
    # src/dynamical_model.cpp
    # src/kinematic_model.cpp
    # # src/inverse_kinematic.cpp
    # src/link.cpp
    # src/robot.cpp
    src/main.cpp
)

set(HEADERS
    # include/dynamical_model.h
    # include/kinematic_model.h
    # include/inverse_kinematic.h
    include/link.h
    include/robot.h
    include/utils.h
)

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )
if(Eigen3_FOUND)
    message(STATUS "Eigen3 found")
else()
    message(FATAL_ERROR "Eigen3 NOT found")
endif()

find_package(jsoncpp REQUIRED)
#find_path(jsoncpp_INCLUDE_DIR json/json.h)
#+include_directories(${jsoncpp_INCLUDE_DIR})
#message(STATUS ${jsoncpp_INCLUDE_DIR})
if(jsoncpp_FOUND)
    message(STATUS "jsoncpp found")
else()
    message(FATAL_ERROR "jsoncpp NOT found")
endif()



add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_link_libraries(${PROJECT_NAME} jsoncpp)


add_library(${library_name} SHARED ${SOURCES} ${HEADERS})
target_link_libraries(${library_name} jsoncpp)

install(TARGETS ${library_name}
    LIBRARY DESTINATION /usr/local/lib
    ARCHIVE DESTINATION /usr/local/lib
)

install(DIRECTORY include
    DESTINATION /usr/local
    FILES_MATCHING PATTERN "*.h"
)