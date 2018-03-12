cmake_minimum_required(VERSION 2.8.3)

project(r_s_planner)

## Add support for C++11, supported in ROS Kinetic and newer
add_definitions(-std=c++11)

set(CMAKE_CXX_FLAGS "-Wall -Wextra ${CMAKE_CXX_FLAGS}")

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
        roscpp
        roslib
        nav_msgs
        )

## System dependencies are found with CMake's conventions
find_package(Boost REQUIRED)


catkin_package(
        #  INCLUDE_DIRS include
        #  LIBRARIES roscpp
        #  CATKIN_DEPENDS test
        DEPENDS Boost
)

include_directories(include)
include_directories(
        ${catkin_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
)

add_executable(${PROJECT_NAME}_node src/path_node.cpp src/reeds_shepp.cpp src/dubins.cpp)
target_link_libraries(${PROJECT_NAME}_node
        ${catkin_LIBRARIES}
        )
## Declare a C++ library
# add_library(${PROJECT_NAME}
#   src/${PROJECT_NAME}/roscpp.cpp
# )

## Add cmake target dependencies of the library
## as an example, code may need to be generated before libraries
## either from message generation or dynamic reconfigure
# add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Declare a C++ executable
## With catkin_make all packages are built within a single CMake context
## The recommended prefix ensures that target names across packages don't collide
# add_executable(${PROJECT_NAME}_node src/roscpp_node.cpp)

## Rename C++ executable without prefix
## The above recommended prefix causes long target names, the following renames the
## target back to the shorter version for ease of user use
## e.g. "rosrun someones_pkg node" instead of "rosrun someones_pkg someones_pkg_node"
# set_target_properties(${PROJECT_NAME}_node PROPERTIES OUTPUT_NAME node PREFIX "")

## Add cmake target dependencies of the executable
## same as for the library above
# add_dependencies(${PROJECT_NAME}_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Specify libraries to link a library or executable target against
# target_link_libraries(${PROJECT_NAME}_node
#   ${catkin_LIBRARIES}
# )