//
// Created by yunfan on 2021/7/16.
//

#include "poly_traj_utils/traj_utils.hpp"
#include "poly_traj_utils/obvp_solver2.hpp"
#include "poly_traj_utils/poly_visual_utils.hpp"
#include "iostream"
#include "poly_traj_utils/scope_timer.hpp"

using namespace std;
ObvpSolver2 solver(1000, 7, 5);
Vec3 points[2];
Vec3 vels[2];
int idx = 0;
ros::Publisher marker_pub, end_point_state_pub, start_point_state_pub;

void wayPoint_callback(const geometry_msgs::PoseStamped::ConstPtr &msg) {
    Mat33 R = MsgUtils::poseToEigenRotation(msg->pose);
    points[idx] = Vec3(msg->pose.position.x, msg->pose.position.y, msg->pose.position.z);
    ROS_INFO("Get target at %lf,%lf,%lf", points[idx].x(), points[idx].y(), points[idx].z());


    if (points[idx].z() < 1e-3)
        points[idx].z() = 2.5;
    Vec3 W(1, 0, 0);
    vels[idx] = R * W;
    if (idx) {
        Piece cur_p2;
        visualization_msgs::Marker list1, list2, acc_1, acc_2;
        geometry_msgs::PoseStamped point;
        point.pose = msg->pose;
        point.pose.position.z = points[idx].z();
        point.header.frame_id = "world";
        point.header.stamp = ros::Time::now();
        end_point_state_pub.publish(point);

        StatePVAJ start_state, end_state;
        start_state << points[0], vels[0], 0, 0, 0, 0, 0, 0;
        end_state << points[1], vels[1], 0, 0, 0, 0, 0, 0;
        Piece pie ;
        {
            TimeConsuming t_("OPT",1e2);
            int cnt = 1e2;
            while (cnt --)
            pie = solver.GenFixPVMinSnapOptT(start_state, end_state);
//            pie = solver.GenFixStateMinSnapOptT(start_state, end_state);
        }
        double vel = pie.getMaxVelRate();
        double acc = pie.getMaxAccRate();
        printf("max_vel = %lf, max_acc = %lf\n", vel, acc);

        double t = solver.GetLastTstar();
        double cost = solver.GetLastCost();
        printf("t = %lf, cost = %lf\n", t, cost);
        PolyVisual::publishPieceToMarker(marker_pub, pie, true);
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

int main(int argc, char **argv) {
    ros::init(argc, argv, "opt_test");
    ros::NodeHandle nh("~");
    ros::Subscriber waypoints_sub = nh.subscribe("/goal", 1, wayPoint_callback);
    ros::Publisher bvp_pub = nh.advertise<nav_msgs::Path>("/path", 1);
    marker_pub = nh.advertise<visualization_msgs::MarkerArray>("/marker", 1);
    start_point_state_pub = nh.advertise<geometry_msgs::PoseStamped>("/start_point_vis", 10);
    end_point_state_pub = nh.advertise<geometry_msgs::PoseStamped>("/end_point_vis", 10);

    ros::AsyncSpinner spinner(0);
    spinner.start();
    ros::Duration(1.0).sleep();

    StatePVAJ start, goal;
    start.head(3) = Vec3(0, 0, 0);
    start.segment(3, 3) = Vec3(1, 2, 1);
    start.segment(6, 3) = Vec3(1, 0, 0);
    start.tail(3) = Vec3(0, 0, 0);

    goal.head(3) = Vec3(5, 0, 0);
    goal.segment(3, 3) = Vec3(2, 0, 1);
    goal.segment(6, 3) = Vec3(0, 1, 0);
    goal.tail(3) = Vec3(0, 0, 0);

//    Piece pie = solver.GenFixStateMinSnapOptTC(start,goal);
////     Time compareson of different
//    {
//        TimeConsuming t__("BVP with t",1e5);
//        int cnt=1e5;
//        while(cnt--)
//        pie = solver.GenFixStateMinSnap(start,goal,5);
//    }
//    {
//        TimeConsuming t__("BVP + Opt t",1e5);
//        int cnt=1e5;
//        while(cnt--)
//
//        pie = solver.GenFixStateMinSnapOptT(start,goal);
//    }
//    {
//        TimeConsuming t__("CBVP",1e5);
//        int cnt=1e5;
//        while(cnt--)
//        pie = solver.GenFixStateMinSnapOptTC(start,goal);
//    }
//    double vel = pie.getMaxVelRate();
//    double acc = pie.getMaxAccRate();
//    printf("max_vel = %lf, max_acc = %lf\n", vel,acc);
//
//    double t = solver.GetLastTstar();
//    double cost = solver.GetLastCost();
//    printf("t = %lf, cost = %lf\n", t,cost);
//    PolyVisual::publishPieceToMarker(marker_pub,pie, true);


//    Piece pie = solver.GenFixStateMinSnapOptT(start, goal);
//    double vel = pie.getMaxVelRate();
//    double acc = pie.getMaxAccRate();
//    printf("max_vel = %lf, max_acc = %lf\n", vel,acc);
//
//    double t = solver.GetLastTstar();
//    double cost = solver.GetLastCost();
//    printf("t = %lf, cost = %lf\n", t,cost);
    Piece pie = solver.GenFixPVMinSnapOptT(start,goal);
    double t = solver.GetLastTstar();
//    double cur_t = 0, acc_cost = 0;
//    double max_v = 0.0;
//    while (cur_t <= t) {
//        Vec3 cur_snap = pie.getSnap(cur_t);
//        acc_cost+=cur_snap.squaredNorm()*0.001;
//        cur_t += 0.001;
//    }
    for(t = 2;t<5;t+=0.01){
        Piece pie = solver.GenFixPVMinSnap(start,goal,t,true);
        fmt::print("{},{};\n",solver.GetLastCost(), t);//3.2550270956988907
    }


    cout<<pie.checkMaxVelRate(2.3)<<endl;



//    double best_t = solver.CalcOptimalDuration(start, goal); //3.1
//    for (double t = best_t - 1; t < best_t + 2; t += 0.05) {
//        Piece pie = solver.GenFixStateMinSnap(start, goal, t);
//        double cur_t = 0, acc_cost = 0;
//        while (cur_t < best_t) {
//            Vec3 cur_snap = pie.getSnap(cur_t);
//            acc_cost += (cur_snap.norm() * cur_snap.norm()) * 0.001;
//            cur_t += 0.001;
//        }
//        double cost = solver.EvaluateSnapCost(pie);
//        printf("%lf, %lf, %lf;\n", t,acc_cost+t*1000,cost);
//        PolyVisual::publishPieceToPath(bvp_pub, pie);
//        if(t - best_t<0.1)
//            PolyVisual::publishPieceToMarker(marker_pub,pie, true);
//    }


    ros::waitForShutdown();
    return 0;
}