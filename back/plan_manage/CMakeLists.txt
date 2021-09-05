cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)
project(plan_manage)

set(CMAKE_BUILD_TYPE "Release")
ADD_COMPILE_OPTIONS(-std=c++14)
ADD_COMPILE_OPTIONS(-std=c++14)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")
#add_compile_options(-fPIC)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# set(CUDA_SEPARABLE_COMPILATION ON)
find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(PCL 1.7 REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem)
find_package(catkin REQUIRED COMPONENTS
        roscpp
        roslib
        rospy
        std_msgs
        visualization_msgs
        cv_bridge
        path_searching
        decomp_ros_utils
        )


catkin_package(
        INCLUDE_DIRS include
        LIBRARIES path_searching
        DEPENDS OpenCV Eigen Boost

        #  DEPENDS system_lib
)
include_directories(
        SYSTEM
        include
        ${catkin_INCLUDE_DIRS}
        ${PROJECT_SOURCE_DIR}/include
        ${Eigen_INCLUDE_DIRS}
        ${Eigen3_INCLUDE_DIRS}
        ${PCL_INCLUDE_DIRS}
        ${DECOMP_UTIL_INCLUDE_DIRS}

)

set(THIRD_PARTY
        ${catkin_LIBRARIES}
        ${PCL_LIBRARIES}
        ${OpenCV_LIBRARIES}
        ${Boost_LIBRARIES}
        -lfmt
        )


add_executable(se3_node
        src/se3_node.cpp
        src/se3_planner.cpp
        )
target_link_libraries(se3_node
        ${THIRD_PARTY}
        )

add_library(plan_manage
        src/se3_planner.cpp
        )
target_link_libraries(plan_manage
        ${THIRD_PARTY}
        )

add_executable(planner_test
        Apps/planner_test.cpp
        src/se3_planner.cpp
        )
target_link_libraries(planner_test
        ${THIRD_PARTY}
        )

# add_library()