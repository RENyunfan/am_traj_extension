//
// Created by yunfan on 2021/3/19.
//

#ifndef SRC_POLY_VISUAL_UTILS_HPP
#define SRC_POLY_VISUAL_UTILS_HPP

#include "ros/ros.h"
#include "visualization_msgs/Marker.h"
#include "visualization_msgs/MarkerArray.h"
#include "geometry_msgs/PoseStamped.h"
#include "nav_msgs/Path.h"
#include "nav_msgs/Odometry.h"
#include "sensor_msgs/PointCloud2.h"
#include "poly_data_structure.hpp"
#include "vector"
#include "Eigen/Dense"
#include "tf/tf.h"
#include "fmt/color.h"

using namespace std;
using namespace Eigen;
typedef Matrix<double, 3, 3> Mat33;
typedef Matrix<double, 3, 1> Vec3;
typedef  std::pair<Vec3,double> Ball;

class PolyVisual{
public:

    static void publishBallToMarkers(ros::Publisher & publisher_, vector<Ball> & searched_pts,int & id){

        ros::Rate loopRate(5000);
        visualization_msgs::Marker  ball_list;
        ball_list.action = visualization_msgs::Marker::DELETEALL;
        publisher_.publish(ball_list);
        for(size_t i = 0 ; i < searched_pts.size(); i++){
            Vec3 cur_pos = searched_pts[i].first;
            double D = searched_pts[i].second;
            ball_list.header.frame_id = "world";
            ball_list.header.stamp = ros::Time::now();
            ball_list.ns = "SFC";
            ball_list.id = id++;
            ball_list.action = visualization_msgs::Marker::ADD;
            ball_list.pose.orientation.w = 1.0;
            ball_list.type = visualization_msgs::Marker::SPHERE;
            // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
            ball_list.scale.x = D*2;
            ball_list.scale.y = D*2;
            ball_list.scale.z = D*2;
            // Line list is red
            ball_list.color.r = 0.5;
            ball_list.color.g = 0.6;
            ball_list.color.b = 1.0;
            ball_list.color.a = 0.6;
            // Create the vertices for the points and lines
            geometry_msgs::Point p;
            p.x = cur_pos.x();
            p.y = cur_pos.y();
            p.z = cur_pos.z();
            // The line list needs two points for each line
            ball_list.pose.position = p;
//                ball_list.points.push_back(p);
            publisher_.publish(ball_list);
            loopRate.sleep();
        }
    }

    static void publishVectorPointsToMarker(ros::Publisher & publisher_, vector<Vec3> & points){
        ros::Rate loopRate(5000);
        visualization_msgs::Marker  ball_list;
        ball_list.action = visualization_msgs::Marker::DELETEALL;
        publisher_.publish(ball_list);
        int id = 0;
        for(size_t i = 0 ; i < points.size(); i++){
            Vec3 cur_pos = points[i];
            ball_list.header.frame_id = "world";
            ball_list.header.stamp = ros::Time::now();
            ball_list.ns = "waypts";
            ball_list.id = id++;
            ball_list.action = visualization_msgs::Marker::ADD;
            ball_list.pose.orientation.w = 1.0;
            ball_list.type = visualization_msgs::Marker::SPHERE;
            // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
            ball_list.scale.x = 0.3;
            ball_list.scale.y = 0.3;
            ball_list.scale.z = 0.3;
            // Line list is red
            ball_list.color.r = 0.5;
            ball_list.color.g = 0.6;
            ball_list.color.b = 1.0;
            ball_list.color.a = 1.0;
            // Create the vertices for the points and lines
            geometry_msgs::Point p;
            p.x = cur_pos.x();
            p.y = cur_pos.y();
            p.z = cur_pos.z();
            // The line list needs two points for each line
            ball_list.pose.position = p;
//                ball_list.points.push_back(p);
            publisher_.publish(ball_list);
            loopRate.sleep();
        }
    }

    static void publishVectorPointsToOdom(ros::Publisher & publisher_, vector<Vec3> & points){
        nav_msgs::Odometry odom;
        ros::Rate loop_rate(1000);
        for( auto i = 0 ; i < points.size(); i++){
            odom.pose.pose.position.x = points[i].x();
            odom.pose.pose.position.y = points[i].y();
            odom.pose.pose.position.z = points[i].z();
            odom.header.frame_id  = "world";
            odom.header.stamp = ros::Time::now();
            publisher_.publish(odom);
            loop_rate.sleep();
        }

    }

    /*--------------  Piece  -------------------------------------*/
    static void publishPieceToPath(ros::Publisher & pub_, Piece & piece){
        if(piece.getDuration() <1e-3 || piece.empty() || piece.getDuration()>999)
            return;
//        fmt::print("Before check the piece\n");
//        fmt::print(fg(fmt::color::azure),"The length of piece size is {}\n",piece.getLength());
//        fmt::print(fg(fmt::color::azure),"The points size is {}\n", points.size());
        vector<Vec3> points = piece.getTraj(0.01);
        nav_msgs::Path pub_path;
        pub_path.poses.reserve(int(points.size()*1.5));
        geometry_msgs::PoseStamped cur_point;
        for( size_t i = 0 ; i < points.size(); i++){
            cur_point.pose.position.x = points[i].x();
            cur_point.pose.position.y = points[i].y();
            cur_point.pose.position.z = points[i].z();
            pub_path.poses.push_back(cur_point);
            pub_path.header.frame_id = "world";
        }
        pub_path.header.frame_id = "world";
        pub_path.header.stamp = ros::Time::now();
        pub_.publish(pub_path);
    }

    static void publishPieceVelToMarkerArray(ros::Publisher &pub_, Piece piece) {
        double t_sum = piece.getDuration();
        double sum_dis = piece.getLength();
        double eval_t = 0.0;
        double v_n;
        visualization_msgs::MarkerArray markers;
        Vec3 cur_vel_dir, cur_vel, cur_pos, end_point;
        Vec3 v2 = (piece.getPos(t_sum - 1e-3) - piece.getPos(0)).normalized();
        int cnt = 0;
        while (eval_t + 1e-4 < t_sum) {
            cur_vel = piece.getVel(eval_t);
            cur_pos = piece.getPos(eval_t);
            v_n = cur_vel.norm();
            cur_vel_dir = cur_vel.cross(
                    Vec3(0, 0, 1)).normalized();// (rot_ang.toRotationMatrix()*cur_vel ).normalized();
                    end_point = v_n * cur_vel_dir;
                    {
                        visualization_msgs::Marker line_list;
                        line_list.header.frame_id = "world";
                        line_list.header.stamp = ros::Time::now();
                        line_list.ns = "vels";
                        line_list.id = cnt++;
                        line_list.action = visualization_msgs::Marker::ADD;
                        line_list.pose.orientation.w = 1.0;
                        line_list.type = visualization_msgs::Marker::LINE_LIST;
                        // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                        line_list.scale.x = 0.03;
                        // Line list is red
                        line_list.color.r = 1.0;
                        line_list.color.a = 1.0;
                        // Create the vertices for the points and lines
                        geometry_msgs::Point p;
                        p.x = cur_pos.x();
                        p.y = cur_pos.y();
                        p.z = cur_pos.z();
                        // The line list needs two points for each line
                        line_list.points.push_back(p);
                        p.x += end_point.x();
                        p.y += end_point.y();
                        p.z += end_point.z();
                        line_list.points.push_back(p);
                        markers.markers.push_back(line_list);
                    }
                    eval_t += t_sum / sum_dis / 8;
        }
        pub_.publish(markers);
    }

    static void publishPieceAccToMarkerArray(ros::Publisher &pub_, Piece piece) {
        double t_sum = piece.getDuration();
        double sum_dis = piece.getLength();
        double eval_t = 0.0;
        double v_n;
        Vec3 cur_vel_dir, cur_acc, cur_pos, end_point, cur_vel;
        Vec3 v2 = (piece.getPos(t_sum - 1e-3) - piece.getPos(0)).normalized();
        visualization_msgs::MarkerArray markers;
        int cnt = 0;
        while (eval_t + 1e-4 < t_sum) {
            cur_acc = piece.getAcc(eval_t);
            cur_vel = piece.getVel(eval_t);
            cur_pos = piece.getPos(eval_t);
            v_n = cur_acc.norm();

            cur_vel_dir = -cur_vel.cross(
                    Vec3(0, 0, 1)).normalized();// (rot_ang.toRotationMatrix()*cur_vel ).normalized();
                    end_point = v_n * cur_vel_dir;
                    {
                        visualization_msgs::Marker line_list;
                        line_list.header.frame_id = "world";
                        line_list.header.stamp = ros::Time::now();
                        line_list.ns = "acc";
                        line_list.id = cnt++;
                        line_list.action = visualization_msgs::Marker::ADD;
                        line_list.pose.orientation.w = 1.0;
                        line_list.type = visualization_msgs::Marker::LINE_LIST;
                        // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                        line_list.scale.x = 0.03;
                        // Line list is blue
                        line_list.color.b = 1.0;
                        line_list.color.a = 1.0;
                        // Create the vertices for the points and lines
                        geometry_msgs::Point p;
                        p.x = cur_pos.x();
                        p.y = cur_pos.y();
                        p.z = cur_pos.z();
                        // The line list needs two points for each line
                        line_list.points.push_back(p);
                        p.x += end_point.x();
                        p.y += end_point.y();
                        p.z += end_point.z();
                        line_list.points.push_back(p);

                        markers.markers.push_back(line_list);
                    }
                    eval_t += t_sum / sum_dis / 8;
        }
        pub_.publish(markers);
    }

    /*--------------  TRAJ  -------------------------------------*/
    static void publishTrajToPath(ros::Publisher & pub_, Trajectory traj){
        if(traj.getTotalDuration() <1e-3)
            return;
        vector<Vec3> points;
        int num_seg = traj.getPieceNum();
        for(size_t i = 0 ; i < num_seg; i++){
            vector<Vec3> cur_points = traj[i].getTraj(0.01);
            points.insert(points.end(),cur_points.begin(),cur_points.end());
        }
        nav_msgs::Path pub_path;
        geometry_msgs::PoseStamped cur_point;
        for( size_t i = 0 ; i < points.size(); i++){
            cur_point.pose.position.x = points[i].x();
            cur_point.pose.position.y = points[i].y();
            cur_point.pose.position.z = points[i].z();
            pub_path.poses.push_back(cur_point);
            pub_path.header.frame_id = "world";
        }
        pub_path.header.frame_id = "world";
        pub_path.header.stamp = ros::Time::now();
        pub_.publish(pub_path);
    }

    static void publishTrajVelToMarker(ros::Publisher & pub_, Trajectory traj, visualization_msgs::Marker & line_list){
        ros::Rate loopRate(5000);
        double t_sum = traj.getTotalDuration();
        double sum_dis = traj.getTotalLength();
        double eval_t = 0.0; double v_n;
        Vec3 cur_vel_dir,cur_vel,cur_pos,end_point;
        Vec3 v2 = (traj.getPos(-1) - traj.getPos(0) ).normalized();
        while(eval_t +1e-4 < t_sum ){
            cur_vel = traj.getVel(eval_t);
            cur_pos = traj.getPos(eval_t);
            v_n = cur_vel.norm();
            cur_vel_dir =cur_vel.cross(Vec3(0,0,1)).normalized();// (rot_ang.toRotationMatrix()*cur_vel ).normalized();
            end_point = v_n * cur_vel_dir;
            {
                line_list.header.frame_id = "world";
                line_list.header.stamp = ros::Time::now();
                line_list.ns = "vels";
                line_list.id = 2;
                line_list.action = visualization_msgs::Marker::ADD;
                line_list.pose.orientation.w = 1.0;
                line_list.type = visualization_msgs::Marker::LINE_LIST;
                // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                line_list.scale.x = 0.03;
                // Line list is red
                line_list.color.r = 1.0;
                line_list.color.a = 1.0;
                // Create the vertices for the points and lines
                geometry_msgs::Point p;
                p.x = cur_pos.x();
                p.y = cur_pos.y();
                p.z = cur_pos.z();
                // The line list needs two points for each line
                line_list.points.push_back(p);
                p.x += end_point.x();
                p.y += end_point.y();
                p.z += end_point.z();
                line_list.points.push_back(p);
                pub_.publish(line_list);
            }
            eval_t += t_sum / sum_dis/8;
            loopRate.sleep();
        }
    }

    static void publishTrajAccToMarker(ros::Publisher & pub_, Trajectory traj, visualization_msgs::Marker & line_list){
        ros::Rate loopRate(5000);
        double t_sum = traj.getTotalDuration();
        double sum_dis = traj.getTotalLength();
        double eval_t = 0.0; double v_n;
        Vec3 cur_vel_dir,cur_acc,cur_pos,end_point,cur_vel;
        Vec3 v2 = (traj.getPos(-1) - traj.getPos(0) ).normalized();

        while(eval_t +1e-4 < t_sum ){
            cur_acc = traj.getAcc(eval_t);
            cur_vel = traj.getVel(eval_t);
            cur_pos = traj.getPos(eval_t);
            v_n = cur_acc.norm();
            cur_vel_dir =-cur_vel.cross(Vec3(0,0,1)).normalized();// (rot_ang.toRotationMatrix()*cur_vel ).normalized();
            end_point = v_n * cur_vel_dir;
            {
                line_list.header.frame_id = "world";
                line_list.header.stamp = ros::Time::now();
                line_list.ns = "acc";
                line_list.id = 2;
                line_list.action = visualization_msgs::Marker::ADD;
                line_list.pose.orientation.w = 1.0;
                line_list.type = visualization_msgs::Marker::LINE_LIST;
                // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                line_list.scale.x = 0.03;
                // Line list is blue
                line_list.color.b = 1.0;
                line_list.color.a = 1.0;
                // Create the vertices for the points and lines
                geometry_msgs::Point p;
                p.x = cur_pos.x();
                p.y = cur_pos.y();
                p.z = cur_pos.z();
                // The line list needs two points for each line
                line_list.points.push_back(p);
                p.x += end_point.x();
                p.y += end_point.y();
                p.z += end_point.z();
                line_list.points.push_back(p);
                pub_.publish(line_list);
            }
            eval_t += t_sum / sum_dis/8;
            loopRate.sleep();
        }
    }

    static void publishTrajJerkToMarker(ros::Publisher & pub_, Trajectory traj, visualization_msgs::Marker & line_list){
        double t_sum = traj.getTotalDuration();
        double sum_dis = traj.getTotalLength();
        double eval_t = 0.0; double v_n;
        Vec3 cur_vel_dir,cur_jerk,cur_pos,end_point,cur_vel;
        ros::Rate loopRate(5000);
        Vec3 v2 = (traj.getPos(-1) - traj.getPos(0) ).normalized();

        while(eval_t +1e-4 < t_sum ){
            cur_jerk = traj.getJerk(eval_t);
            cur_vel = traj.getVel(eval_t);
            cur_pos = traj.getPos(eval_t);
            v_n = cur_jerk.norm();
            cur_vel_dir =-cur_vel.cross(Vec3(0,0,1)).normalized();// (rot_ang.toRotationMatrix()*cur_vel ).normalized();
            end_point = v_n * cur_vel_dir;
            {
                line_list.header.frame_id = "world";
                line_list.header.stamp = ros::Time::now();
                line_list.ns = "jerk";
                line_list.id = 2;
                line_list.action = visualization_msgs::Marker::ADD;
                line_list.pose.orientation.w = 1.0;
                line_list.type = visualization_msgs::Marker::LINE_LIST;
                // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                line_list.scale.x = 0.03;
                // Line list is blue
                line_list.color.g = 1.0;
                line_list.color.a = 1.0;
                // Create the vertices for the points and lines
                geometry_msgs::Point p;
                p.x = cur_pos.x();
                p.y = cur_pos.y();
                p.z = cur_pos.z();
                // The line list needs two points for each line
                line_list.points.push_back(p);
                p.x += end_point.x();
                p.y += end_point.y();
                p.z += end_point.z();
                line_list.points.push_back(p);
                pub_.publish(line_list);
            }
            eval_t += t_sum / sum_dis/8;
            loopRate.sleep();
        }
    }

};

class MsgUtils {
public:
    static Mat33 poseToEigenRotation(geometry_msgs::Pose msg) {
        Mat33 R_eig;
        tf::Quaternion quat;
        tf::quaternionMsgToTF(msg.orientation, quat);
        tf::Matrix3x3 R = tf::Matrix3x3(quat);
        R_eig << R[0][0], R[0][1], R[0][2],
                R[1][0], R[1][1], R[1][2],
                R[2][0], R[2][1], R[2][2];
        return R_eig;
    }

    static Vec3 poseToEigenVec(geometry_msgs::Pose pose) {
        return Vec3(pose.position.x, pose.position.y, pose.position.z);
    }

    static double deg2rad(double deg){
        return deg*3.141592/180.0;
    }
};


#endif //SRC_POLY_VISUAL_UTILS_HPP