cmake_minimum_required(VERSION 3.10)

project(test)

find_package(orocos_kdl REQUIRED)

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )

find_package(catkin REQUIRED COMPONENTS roscpp)
include_directories(${catkin_INCLUDE_DIRS})

find_library(URDF_LIBRARY NAMES urdf PATHS /usr/local/lib/liburdf.so)
if(URDF_LIBRARY)
    
else()
    message(FATAL_ERROR "liburdf not found")
endif()

find_library(KDL_PARSER_LIBRARY NAMES kdl_parser PATHS /usr/lib/libkdl_parser.so)
if(KDL_PARSER_LIBRARY)
    
else()
    message(FATAL_ERROR "libkdl_parser not found")
endif()


include_directories(
  ${catkin_INCLUDE_DIRS}
)
list(APPEND CMAKE_PREFIX_PATH "/opt/ros/melodic/lib/")
find_library(ROS_CONSOLE_BRIDGE_LIB NAMES rosconsole_bridge PATHS /opt/ros/melodic/lib/)
if(ROS_CONSOLE_BRIDGE_LIB)
    message(STATUS "Libreria rosconsole_bridge trovata: ${ROS_CONSOLE_BRIDGE_LIB}")
else()
    message(FATAL_ERROR "Libreria rosconsole_bridge non trovata")
endif()


add_executable(${PROJECT_NAME} dynamics.cpp)
target_link_libraries(${PROJECT_NAME} ${orocos_kdl_LIBRARIES} ${catkin_LIBRARIES} ${URDF_LIBRARY} ${KDL_PARSER_LIBRARY} ${ROS_CONSOLE_BRIDGE_LIB})