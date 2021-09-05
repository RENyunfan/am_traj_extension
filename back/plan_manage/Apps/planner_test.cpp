//
// Created by yunfan on 2021/8/18.
//

#include <plan_manage/se3_planner.h>
#include "scope_timer.hpp"
#include "msg_utils.h"

MavGlobalPlanner::Ptr planner_ptr;
ros::Publisher traj_snap_pub, traj_opt_pub, end_point_state_pub, start_point_state_pub, sp_pt_pub;
Vec3 points[2];
Vec3 vels[2];
Vec3 W(2, 0, 0);
int idx = 0;

void WaypointCallback(const geometry_msgs::PoseStamped::ConstPtr &msg) {

    Mat33 R = MsgUtils::poseToEigenRotation(msg->pose);
    points[idx] = Vec3(msg->pose.position.x, msg->pose.position.y, msg->pose.position.z);
    ROS_INFO("Get target at %lf,%lf,%lf", points[idx].x(), points[idx].y(), points[idx].z());

    points[idx].z() = 1;
    Vec3 W(3, 0, 0);
    vels[idx] = R * W;
    if (idx) {
        geometry_msgs::PoseStamped point;
        point.pose = msg->pose;
        point.pose.position.z = points[idx].z();
        point.header.frame_id = "world";
        point.header.stamp = ros::Time::now();
        end_point_state_pub.publish(point);
        Mat33 start, goal;
        start << points[0], vels[0], Vec3(0, 0, 0);
        goal << points[1], vels[1], Vec3(0, 0, 0);

        planner_ptr->plan(start, goal);

        idx = 0;
    } else {
        geometry_msgs::PoseStamped point;
        point.pose = msg->pose;
        point.pose.position.z = points[idx].z();
        point.header.frame_id = "world";
        point.header.stamp = ros::Time::now();
        start_point_state_pub.publish(point);
        idx = 1;
    }

}

void OdomCallback(const geometry_msgs::TwistStampedConstPtr & msg){
    Vec3 vel = Vec3(msg->twist.linear.x,msg->twist.linear.y,msg->twist.linear.z);
    fmt::print("{}\n",vel.norm());
}

int main(int argc, char **argv) {
    ros::init(argc, argv, "planner_test");
    ros::NodeHandle nh("~");
    Config config;
    config.loadParameters(nh);
    planner_ptr.reset(new MavGlobalPlanner(config, nh));
    start_point_state_pub = nh.advertise<geometry_msgs::PoseStamped>("/start_pt", 1);
    end_point_state_pub = nh.advertise<geometry_msgs::PoseStamped>("/end_pt", 1);
    ros::Subscriber goal_sub = nh.subscribe("/goal", 1, WaypointCallback);
    ros::Subscriber odom_sub = nh.subscribe("/lio/vision_vel/twist",5,OdomCallback);
    ros::AsyncSpinner spinner(0);
    spinner.start();
    ros::Duration(0.1).sleep();


    ros::waitForShutdown();
    return 0;
}