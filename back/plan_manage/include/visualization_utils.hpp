//
// Created by yunfan on 2021/6/4.
//
/*
 * SRC_VISUALIZATION_UTILS_H
 *
 *
 * */

#ifndef SRC_VISUALIZATION_UTILS_H
#define SRC_VISUALIZATION_UTILS_H

#include "ros/ros.h"
#include "visualization_msgs/Marker.h"
#include "visualization_msgs/MarkerArray.h"
#include "geometry_msgs/PoseStamped.h"
#include "nav_msgs/Path.h"
#include "nav_msgs/Odometry.h"
#include "sensor_msgs/PointCloud2.h"
#include "vector"
#include "Eigen/Dense"
#include "tf/tf.h"
#include "cstring"
#include "se3gcopter/se3gcopter.hpp"

using namespace std;
using namespace Eigen;
typedef Matrix<double, 3, 3> Mat33;
typedef Matrix<double, 3, 1> Vec3;

#define RED (Vec3(255,0,0))
#define CORNBLUE (Vec3(100,149,237))
#define BLUE (Vec3(0,0,237))
#define GREENYELLOW (Vec3(173,255,47))
#define SPRINGGREEN (Vec3(0 ,139 ,69))
#define INDEANRED (Vec3(255, 106, 106))

#define ORANGE (Vec3(255, 165, 0))
#define YELLOW (Vec3(255, 215, 0))
#define RBLUE (    Vec3(39 ,64 ,139))

class VisualUtils {
private:
    static ros::Publisher nullpub;
public:
    /* static void PublishVectorVec3ToMarkerArray(ros::Publisher &publisher_, vector<Vec3> &points,
                                                Vec3 color = Vec3(0.5, 0.5, 0.5), string pt_ns = "vecpts") {
         visualization_msgs::MarkerArray mrkarr;
         color /= 255;
         visualization_msgs::Marker mk;
         mk.action = visualization_msgs::Marker::DELETEALL;
         mrkarr.markers.push_back(mk);

         int id = 0;
         for (size_t i = 0; i < points.size(); i++) {
             visualization_msgs::Marker ball_list;
             Vec3 cur_pos = points[i];
             ball_list.header.frame_id = "world";
             ball_list.header.stamp = ros::Time::now();
             ball_list.ns = pt_ns.c_str();
             ball_list.id = id++;
             ball_list.action = visualization_msgs::Marker::ADD;
             ball_list.pose.orientation.w = 1.0;
             ball_list.type = visualization_msgs::Marker::SPHERE;
             // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
             ball_list.scale.x = 0.3;
             ball_list.scale.y = 0.3;
             ball_list.scale.z = 0.3;
             // Line list is red
             ball_list.color.r = color[0];
             ball_list.color.g = color[1];
             ball_list.color.b = color[2];
             ball_list.color.a = 1.0;
             // Create the vertices for the points and lines
             geometry_msgs::Point p;
             p.x = cur_pos.x();
             p.y = cur_pos.y();
             p.z = cur_pos.z();
             // The line list needs two points for each line
             ball_list.pose.position = p;
             //                ball_list.points.push_back(p);
             mrkarr.markers.push_back(ball_list);
         }
         publisher_.publish(mrkarr);
     }

     static void PublishPieceToMarkerArray(ros::Publisher &pub_, Piece piece, string ns = "line", Vec3 color = Vec3(1, 0, 0)) {
         visualization_msgs::MarkerArray mrkarr;
         visualization_msgs::Marker mk;
         color /= 255;
         static int line_cnt = 0;

         if (piece.getDuration() < 1e-3) {
             return;
         }
         vector<Vec3> points = piece.getTraj(0.01);
         nav_msgs::Path pub_path;
         geometry_msgs::PoseStamped cur_point;
         for (size_t i = 0; i < points.size() - 1; i++) {
             geometry_msgs::Point p;
             // publish lines
             visualization_msgs::Marker line_list;
             line_list.header.frame_id = "world";
             line_list.header.stamp = ros::Time::now();
             line_list.ns = ns.c_str();
             line_list.id = line_cnt++;
             line_list.action = visualization_msgs::Marker::ADD;
             line_list.pose.orientation.w = 1.0;
             line_list.type = visualization_msgs::Marker::LINE_LIST;
             // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
             line_list.scale.x = 0.03;
             // Line list is blue
             line_list.color.r = color[0];
             line_list.color.g = color[1];
             line_list.color.b = color[2];
             line_list.color.a = 1.0;
             // Create the vertices for the points and lines

             p.x = points[i].x();
             p.y = points[i].y();
             p.z = points[i].z();
             // The line list needs two points for each line
             line_list.points.push_back(p);
             p.x = points[i + 1].x();
             p.y = points[i + 1].y();
             p.z = points[i + 1].z();
             // The line list needs
             line_list.points.push_back(p);

             mrkarr.markers.push_back(line_list);
         }
         pub_.publish(mrkarr);

     }

     static void PublishPieceFuallStateToMarkerArray(ros::Publisher &pub_, Piece piece) {

         VisualUtils::PublishPiecePosiToMarkerArray(pub_, piece);
         VisualUtils::PublishPieceVelToMarkerArray(pub_, piece);
         VisualUtils::PublishPieceAccToMarkerArray(pub_, piece);
         VisualUtils::PublishPieceJerkToMarkerArray(pub_, piece);
         VisualUtils::PublishPieceSnapToMarkerArray(pub_, piece);
     }

     static void PublishPiecePosiToMarkerArray(ros::Publisher &pub_, Piece piece, Vec3 color = RBLUE) {
         double t_sum = piece.getDuration();
         double sum_dis = piece.getLength();
         double eval_t = 0.01;
         double v_n;
         color /= 255;
         Vec3 cur_vel_dir, cur_acc, cur_pos, end_point, cur_vel;
         Vec3 v2 = (piece.getPos(t_sum - 1e-3) - piece.getPos(0)).normalized();
         visualization_msgs::MarkerArray markers;
         int cnt = 0;
         Vec3 last_posi;
         cur_pos = piece.getPos(0);
         while (eval_t + 1e-4 < t_sum) {
             last_posi = cur_pos;
             cur_pos = piece.getPos(eval_t);
             visualization_msgs::Marker line_list;
             line_list.header.frame_id = "world";
             line_list.header.stamp = ros::Time::now();
             line_list.ns = "posi";
             line_list.id = cnt++;
             line_list.action = visualization_msgs::Marker::ADD;
             line_list.pose.orientation.w = 1.0;
             line_list.type = visualization_msgs::Marker::LINE_LIST;
             // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
             line_list.scale.x = 0.03;
             // Line list is blue
             line_list.color.r = color[0];
             line_list.color.g = color[1];
             line_list.color.b = color[2];
             line_list.color.a = 1;
             // Create the vertices for the points and lines
             geometry_msgs::Point p;
             p.x = last_posi.x();
             p.y = last_posi.y();
             p.z = last_posi.z();
             // The line list needs two points for each line
             line_list.points.push_back(p);
             p.x = cur_pos.x();
             p.y = cur_pos.y();
             p.z = cur_pos.z();
             line_list.points.push_back(p);

             markers.markers.push_back(line_list);
             eval_t += t_sum / sum_dis / 8;
         }
         pub_.publish(markers);
     }

     static void PublishPieceJerkToMarkerArray(ros::Publisher &pub_, Piece piece, Vec3 color = YELLOW) {
         double t_sum = piece.getDuration();
         double sum_dis = piece.getLength();
         double eval_t = 0.0;
         double v_n;
         color /= 255;
         Vec3 cur_vel_dir, cur_acc, cur_pos, end_point, cur_vel;
         Vec3 v2 = (piece.getPos(t_sum - 1e-3) - piece.getPos(0)).normalized();
         visualization_msgs::MarkerArray markers;
         int cnt = 0;
         while (eval_t + 1e-4 < t_sum) {
             cur_acc = piece.getJerk(eval_t);
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
                 line_list.ns = "jerk";
                 line_list.id = cnt++;
                 line_list.action = visualization_msgs::Marker::ADD;
                 line_list.pose.orientation.w = 1.0;
                 line_list.type = visualization_msgs::Marker::LINE_LIST;
                 // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                 line_list.scale.x = 0.03;
                 // Line list is blue
                 line_list.color.r = color[0];
                 line_list.color.g = color[1];
                 line_list.color.b = color[2];
                 line_list.color.a = 0.8;
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

     static void PublishPieceSnapToMarkerArray(ros::Publisher &pub_, Piece piece, Vec3 color = BLUE) {
         double t_sum = piece.getDuration();
         double sum_dis = piece.getLength();
         double eval_t = 0.0;
         double v_n;
         color /= 255;
         Vec3 cur_vel_dir, cur_acc, cur_pos, end_point, cur_vel;
         Vec3 v2 = (piece.getPos(t_sum - 1e-3) - piece.getPos(0)).normalized();
         visualization_msgs::MarkerArray markers;
         int cnt = 0;
         while (eval_t + 1e-4 < t_sum) {
             cur_acc = piece.getSnap(eval_t);
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
                 line_list.ns = "snap";
                 line_list.id = cnt++;
                 line_list.action = visualization_msgs::Marker::ADD;
                 line_list.pose.orientation.w = 1.0;
                 line_list.type = visualization_msgs::Marker::LINE_LIST;
                 // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                 line_list.scale.x = 0.03;
                 // Line list is blue
                 line_list.color.r = color[0];
                 line_list.color.g = color[1];
                 line_list.color.b = color[2];
                 line_list.color.a = 0.8;
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

     static void PublishPieceVelToMarkerArray(ros::Publisher &pub_, Piece piece, Vec3 color = RED) {
         double t_sum = piece.getDuration();
         double sum_dis = piece.getLength();
         double eval_t = 0.0;
         double v_n;
         color /= 255;
         Vec3 cur_vel_dir, cur_acc, cur_pos, end_point, cur_vel;
         Vec3 v2 = (piece.getPos(t_sum - 1e-3) - piece.getPos(0)).normalized();
         visualization_msgs::MarkerArray markers;
         int cnt = 0;
         while (eval_t + 1e-4 < t_sum) {
             cur_acc = piece.getVel(eval_t);
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
                 line_list.ns = "vel";
                 line_list.id = cnt++;
                 line_list.action = visualization_msgs::Marker::ADD;
                 line_list.pose.orientation.w = 1.0;
                 line_list.type = visualization_msgs::Marker::LINE_LIST;
                 // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                 line_list.scale.x = 0.03;
                 // Line list is blue
                 line_list.color.r = color[0];
                 line_list.color.g = color[1];
                 line_list.color.b = color[2];
                 line_list.color.a = 0.8;
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

     static void PublishPieceAccToMarkerArray(ros::Publisher &pub_, Piece piece, Vec3 color = ORANGE) {
         double t_sum = piece.getDuration();
         double sum_dis = piece.getLength();
         double eval_t = 0.0;
         double v_n;
         color /= 255;
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
                 line_list.color.r = color[0];
                 line_list.color.g = color[1];
                 line_list.color.b = color[2];
                 line_list.color.a = 0.8;
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


     static void PublishPieceToPath(ros::Publisher &pub_, Piece piece, string topic = "world") {
         if (piece.getDuration() < 1e-3)
             return;
         vector<Vec3> points = piece.getTraj(0.01);
         nav_msgs::Path pub_path;
         geometry_msgs::PoseStamped cur_point;
         for (size_t i = 0; i < points.size(); i++) {
             cur_point.pose.position.x = points[i].x();
             cur_point.pose.position.y = points[i].y();
             cur_point.pose.position.z = points[i].z();
             pub_path.poses.push_back(cur_point);
             pub_path.header.frame_id = topic.c_str();
         }
         pub_path.header.frame_id = topic.c_str();
         pub_path.header.stamp = ros::Time::now();
         pub_.publish(pub_path);
     }

     static inline void PublishVec3ToMarkerPoint(ros::Publisher publisher_, Vec3 pt, Vec3 color = Vec3(255, 0, 255),
                                                 string name = "test_point") {
         visualization_msgs::Marker marker_ball;
         static int cnt = 0;
         color /= 255;
         marker_ball.action = visualization_msgs::Marker::DELETEALL;
         Vec3 cur_pos = pt;
         marker_ball.header.frame_id = "world";
         marker_ball.header.stamp = ros::Time::now();
         marker_ball.ns = name.c_str();
         marker_ball.id = cnt++;
         marker_ball.action = visualization_msgs::Marker::ADD;
         marker_ball.pose.orientation.w = 1.0;
         marker_ball.type = visualization_msgs::Marker::SPHERE;
         // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
         marker_ball.scale.x = 0.1;
         marker_ball.scale.y = 0.1;
         marker_ball.scale.z = 0.1;
         // Line list is red
         marker_ball.color.r = color[0];
         marker_ball.color.g = color[1];
         marker_ball.color.b = color[2];
         marker_ball.color.a = 0.6;
         // Create the vertices for the points and lines
         geometry_msgs::Point p;
         p.x = cur_pos.x();
         p.y = cur_pos.y();
         p.z = cur_pos.z();
         // The line list needs two points for each line
         marker_ball.pose.position = p;
 //                ball_list.points.push_back(p);
         publisher_.publish(marker_ball);
     }
 */
/*--------------  TRAJ  -------------------------------------*/
    /* static void publishTrajToPath(ros::Publisher &pub_, Trajectory traj) {
         if (traj.getTotalDuration() < 1e-3)
             return;
         vector<Vec3> points;
         int num_seg = traj.getPieceNum();
         for (size_t i = 0; i < num_seg; i++) {
             vector<Vec3> cur_points = traj[i].getTraj(0.01);
             points.insert(points.end(), cur_points.begin(), cur_points.end());
         }
         nav_msgs::Path pub_path;
         geometry_msgs::PoseStamped cur_point;
         for (size_t i = 0; i < points.size(); i++) {
             cur_point.pose.position.x = points[i].x();
             cur_point.pose.position.y = points[i].y();
             cur_point.pose.position.z = points[i].z();
             pub_path.header.frame_id = "world";
             pub_path.poses.push_back(cur_point);
         }
         pub_path.header.frame_id = "world";
         pub_path.header.stamp = ros::Time::now();
         pub_.publish(pub_path);
     }

     static void publishTrajVelToMarkerArray(ros::Publisher &pub_, Trajectory traj, Vec3 color = RED) {
         double t_sum = traj.getTotalDuration();
         double sum_dis = traj.getTotalLength();
         double eval_t = 0.0;
         double v_n;
         Vec3 cur_vel_dir, cur_vel, cur_pos, end_point;
         Vec3 v2 = (traj.getPos(-1) - traj.getPos(0)).normalized();
         visualization_msgs::Marker line_list;
         line_list.action = visualization_msgs::Marker::DELETEALL;
         visualization_msgs::MarkerArray mrkarr;
         mrkarr.markers.push_back(line_list);
         pub_.publish(mrkarr);
         color /= 255;
         mrkarr.markers.clear();
         int idx = 0;
         while (eval_t + 1e-4 < t_sum) {
             visualization_msgs::Marker line_list;
             cur_vel = traj.getVel(eval_t);
             cur_pos = traj.getPos(eval_t);
             v_n = cur_vel.norm();
             printf("%lf;\n", v_n);
             cur_vel_dir = cur_vel.cross(
                     Vec3(0, 0, 1)).normalized();// (rot_ang.toRotationMatrix()*cur_vel ).normalized();
             end_point = v_n * cur_vel_dir;
             {
                 line_list.header.frame_id = "world";
                 line_list.header.stamp = ros::Time::now();
                 line_list.ns = "vels";
                 line_list.id = idx++;
                 line_list.action = visualization_msgs::Marker::ADD;
                 line_list.pose.orientation.w = 1.0;
                 line_list.type = visualization_msgs::Marker::LINE_LIST;
                 // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                 line_list.scale.x = 0.03;
                 // Line list is red
                 line_list.color.r = color[0];
                 line_list.color.g = color[1];
                 line_list.color.b = color[2];
                 line_list.color.a = 0.3;
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

                 mrkarr.markers.push_back(line_list);
             }
             eval_t += t_sum / sum_dis / 32;

         }
         pub_.publish(mrkarr);
     }

     static void publishTrajAccToMarkerArray(ros::Publisher &pub_, Trajectory traj, Vec3 color = BLUE) {
         double t_sum = traj.getTotalDuration();
         double sum_dis = traj.getTotalLength();
         double eval_t = 0.0;
         double a_n;
         color /= 255;
         Vec3 cur_vel_dir, cur_vel, cur_acc, cur_pos, end_point;
         Vec3 v2 = (traj.getPos(-1) - traj.getPos(0)).normalized();
         visualization_msgs::Marker line_list;
         line_list.action = visualization_msgs::Marker::DELETEALL;
         visualization_msgs::MarkerArray mrkarr;
         mrkarr.markers.push_back(line_list);
         pub_.publish(mrkarr);
         mrkarr.markers.clear();
         int idx = 0;
         while (eval_t + 1e-4 < t_sum) {
             visualization_msgs::Marker line_list;
             cur_acc = traj.getAcc(eval_t);
             cur_vel = traj.getVel(eval_t);
             cur_pos = traj.getPos(eval_t);
             a_n = cur_acc.norm();
             cur_vel_dir = -cur_vel.cross(
                     Vec3(0, 0, 1)).normalized();// (rot_ang.toRotationMatrix()*cur_vel ).normalized();
             end_point = a_n * cur_vel_dir + cur_pos;
             {
                 line_list.header.frame_id = "world";
                 line_list.header.stamp = ros::Time::now();
                 line_list.ns = "accs";
                 line_list.id = idx++;
                 line_list.action = visualization_msgs::Marker::ADD;
                 line_list.pose.orientation.w = 1.0;
                 line_list.type = visualization_msgs::Marker::LINE_LIST;
                 // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                 line_list.scale.x = 0.03;
                 // Line list is red
                 line_list.color.r = color[0];
                 line_list.color.g = color[1];
                 line_list.color.b = color[2];
                 line_list.color.a = 0.3;
                 // Create the vertices for the points and lines
                 geometry_msgs::Point p;
                 p.x = cur_pos.x();
                 p.y = cur_pos.y();
                 p.z = cur_pos.z();
                 // The line list needs two points for each line
                 line_list.points.push_back(p);
                 p.x = end_point.x();
                 p.y = end_point.y();
                 p.z = end_point.z();
                 line_list.points.push_back(p);

                 mrkarr.markers.push_back(line_list);
             }
             eval_t += t_sum / sum_dis / 32;

         }
         pub_.publish(mrkarr);
     }
 */
    static void PublishTrajVelToMarkerArray(ros::Publisher &pub_, Trajectory traj, Vec3 color = RED) {
        double t_sum = traj.getTotalDuration();
        double eval_t = 0.0;
        double a_n;
        double sum_dis = (traj.getPos(0) - traj.getPos(t_sum)).norm();
        color /= 255;
        Vec3 cur_vel_dir, cur_vel, cur_acc, cur_pos, end_point;
        Vec3 v2 = (traj.getPos(-1) - traj.getPos(0)).normalized();
        visualization_msgs::Marker line_list;
        line_list.action = visualization_msgs::Marker::DELETEALL;
        visualization_msgs::MarkerArray mrkarr;
        mrkarr.markers.push_back(line_list);
        pub_.publish(mrkarr);
        mrkarr.markers.clear();
        int idx = 0;
        while (eval_t + 1e-4 < t_sum) {
            visualization_msgs::Marker line_list;
            cur_acc = traj.getVel(eval_t);
            cur_vel = traj.getVel(eval_t);
            cur_pos = traj.getPos(eval_t);
            a_n = cur_acc.norm();
            cur_vel_dir = -cur_vel.cross(
                    Vec3(0, 0, 1)).normalized();// (rot_ang.toRotationMatrix()*cur_vel ).normalized();
            end_point = a_n * cur_vel_dir + cur_pos;
            {
                line_list.header.frame_id = "world";
                line_list.header.stamp = ros::Time::now();
                line_list.ns = "vel";
                line_list.id = idx++;
                line_list.action = visualization_msgs::Marker::ADD;
                line_list.pose.orientation.w = 1.0;
                line_list.type = visualization_msgs::Marker::LINE_LIST;
                // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                line_list.scale.x = 0.03;
                // Line list is red
                line_list.color.r = color[0];
                line_list.color.g = color[1];
                line_list.color.b = color[2];
                line_list.color.a = 0.3;
                // Create the vertices for the points and lines
                geometry_msgs::Point p;
                p.x = cur_pos.x();
                p.y = cur_pos.y();
                p.z = cur_pos.z();
                // The line list needs two points for each line
                line_list.points.push_back(p);
                p.x = end_point.x();
                p.y = end_point.y();
                p.z = end_point.z();
                line_list.points.push_back(p);

                mrkarr.markers.push_back(line_list);
            }
            eval_t += t_sum / sum_dis / 32;

        }
        pub_.publish(mrkarr);
    }

    static void PublishTrajJerkToMarkerArray(ros::Publisher &pub_, Trajectory traj, Vec3 color = BLUE) {
        double t_sum = traj.getTotalDuration();
        double eval_t = 0.0;
        double a_n;
        double sum_dis = (traj.getPos(0) - traj.getPos(t_sum)).norm();
        color /= 255;
        Vec3 cur_vel_dir, cur_vel, cur_acc, cur_pos, end_point;
        Vec3 v2 = (traj.getPos(-1) - traj.getPos(0)).normalized();
        visualization_msgs::Marker line_list;
        line_list.action = visualization_msgs::Marker::DELETEALL;
        visualization_msgs::MarkerArray mrkarr;
        mrkarr.markers.push_back(line_list);
        pub_.publish(mrkarr);
        mrkarr.markers.clear();
        int idx = 0;
        while (eval_t + 1e-4 < t_sum) {
            visualization_msgs::Marker line_list;
            cur_acc = traj.getJer(eval_t);
            cur_vel = traj.getVel(eval_t);
            cur_pos = traj.getPos(eval_t);

            a_n = cur_acc.norm();
            cur_vel_dir = -cur_vel.cross(
                    Vec3(0, 0, 1)).normalized();// (rot_ang.toRotationMatrix()*cur_vel ).normalized();
            end_point = a_n * cur_vel_dir + cur_pos;
            {
                line_list.header.frame_id = "world";
                line_list.header.stamp = ros::Time::now();
                line_list.ns = "jerk";
                line_list.id = idx++;
                line_list.action = visualization_msgs::Marker::ADD;
                line_list.pose.orientation.w = 1.0;
                line_list.type = visualization_msgs::Marker::LINE_LIST;
                // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                line_list.scale.x = 0.03;
                // Line list is red
                line_list.color.r = color[0];
                line_list.color.g = color[1];
                line_list.color.b = color[2];
                line_list.color.a = 0.3;
                // Create the vertices for the points and lines
                geometry_msgs::Point p;
                p.x = cur_pos.x();
                p.y = cur_pos.y();
                p.z = cur_pos.z();
                // The line list needs two points for each line
                line_list.points.push_back(p);
                p.x = end_point.x();
                p.y = end_point.y();
                p.z = end_point.z();
                line_list.points.push_back(p);

                mrkarr.markers.push_back(line_list);
            }
            eval_t += t_sum / sum_dis / 32;

        }
        pub_.publish(mrkarr);
    }

    /*   static void PublishTrajSnapToMarkerArray(ros::Publisher &pub_, Trajectory traj, Vec3 color = GREENYELLOW) {
           double t_sum = traj.getTotalDuration();
           double sum_dis = traj.getTotalLength();
           double eval_t = 0.0;
           double a_n;
           color /= 255;
           Vec3 cur_vel_dir, cur_vel, cur_acc, cur_pos, end_point;
           Vec3 v2 = (traj.getPos(-1) - traj.getPos(0)).normalized();
           visualization_msgs::Marker line_list;
           line_list.action = visualization_msgs::Marker::DELETEALL;
           visualization_msgs::MarkerArray mrkarr;
           mrkarr.markers.push_back(line_list);
           pub_.publish(mrkarr);
           mrkarr.markers.clear();
           int idx = 0;
           while (eval_t + 1e-4 < t_sum) {
               visualization_msgs::Marker line_list;
               cur_acc = traj.getSnap(eval_t);
               cur_vel = traj.getVel(eval_t);
               cur_pos = traj.getPos(eval_t);
               a_n = cur_acc.norm();
               cur_vel_dir = -cur_vel.cross(
                       Vec3(0, 0, 1)).normalized();// (rot_ang.toRotationMatrix()*cur_vel ).normalized();
               end_point = a_n * cur_vel_dir + cur_pos;
               {
                   line_list.header.frame_id = "world";
                   line_list.header.stamp = ros::Time::now();
                   line_list.ns = "snap";
                   line_list.id = idx++;
                   line_list.action = visualization_msgs::Marker::ADD;
                   line_list.pose.orientation.w = 1.0;
                   line_list.type = visualization_msgs::Marker::LINE_LIST;
                   // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                   line_list.scale.x = 0.03;
                   // Line list is red
                   line_list.color.r = color[0];
                   line_list.color.g = color[1];
                   line_list.color.b = color[2];
                   line_list.color.a = 0.3;
                   // Create the vertices for the points and lines
                   geometry_msgs::Point p;
                   p.x = cur_pos.x();
                   p.y = cur_pos.y();
                   p.z = cur_pos.z();
                   // The line list needs two points for each line
                   line_list.points.push_back(p);
                   p.x = end_point.x();
                   p.y = end_point.y();
                   p.z = end_point.z();
                   line_list.points.push_back(p);

                   mrkarr.markers.push_back(line_list);
               }
               eval_t += t_sum / sum_dis / 32;

           }
           pub_.publish(mrkarr);
       }*/

};


#endif //SRC_KRRT_PLANNER_H
