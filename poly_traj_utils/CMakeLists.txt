cmake_minimum_required(VERSION 2.8.3)
project(poly_traj_utils)

add_compile_options(-std=c++11)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g -fPIC")

find_package(catkin REQUIRED COMPONENTS
        roscpp
        rospy
        std_msgs
        )

find_package(Eigen3 REQUIRED)

catkin_package(
        INCLUDE_DIRS include
        LIBRARIES poly_traj_utils
        CATKIN_DEPENDS roscpp rospy std_msgs
        #  DEPENDS system_lib
)

include_directories(
        SYSTEM
        include
        ${catkin_INCLUDE_DIRS}
        ${Eigen3_INCLUDE_DIRS}
        ${PCL_INCLUDE_DIRS}
        ${OpenCV_INCLUDE_DIRS}
)

set(THIRD_PARTY
        ${catkin_LIBRARIES}
        ${PCL_LIBRARIES}
        ${OpenCV_LIBRARIES}
        -lfmt
        dw
        )

add_library(${PROJECT_NAME}
        src/poly_traj_utils.cpp
        )
target_link_libraries(
        ${PROJECT_NAME}
        ${THIRD_PARTY}
)
#==========================================================
#
#add_executable(simple_traj_coef
#        Apps/simple_traj_generator_coef.cpp
#        )
#
#target_link_libraries(simple_traj_coef
#        ${PROJECT_NAME}
#        ${THIRD_PARTY}
#        )

#add_executable(poly_utils_test
#        Apps/poly_utils_test.cpp
#        )
#
#target_link_libraries(poly_utils_test
#        ${PROJECT_NAME}
#        ${THIRD_PARTY}
#        )
#
#add_executable(opt_utils_test
#        Apps/opt_utils_test.cpp
#        )
#
#target_link_libraries(opt_utils_test
#        ${PROJECT_NAME}
#        ${THIRD_PARTY}
#        )
add_executable(minco_test
        Apps/minco_test.cpp
        )

target_link_libraries(minco_test
        ${THIRD_PARTY}
        )
