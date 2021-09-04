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
#include "poly_traj_utils/poly_data_structure.hpp"
#include "vector"
#include "Eigen/Dense"
#include "tf/tf.h"
#include "cstring"
#include "nanoflann/state_kdtree.hpp"

using namespace std;
using namespace Eigen;
typedef Matrix<double, 3, 3> Mat33;
typedef Matrix<double, 3, 1> Vec3;
typedef std::pair<Vec3,double> Ball;
typedef vector<std::pair<Vec3,double>> Balls;
#define RED (Vec3(255,0,0))
#define CORNBLUE (Vec3(100,149,237))
#define BLUE (Vec3(0,0,237))
#define GREENYELLOW (Vec3(173,255,47))
#define SPRINGGREEN (Vec3(0 ,139 ,69))
#define INDEANRED (Vec3(255, 106, 106))

#define ORANGE (Vec3(255, 165, 0))
#define YELLOW (Vec3(255, 215, 0))
#define RBLUE (Vec3(39 ,64 ,139))

class Color : public std_msgs::ColorRGBA {
public:
    Color() : std_msgs::ColorRGBA() {}

    Color(double red, double green, double blue) : Color(red, green, blue, 1.0) {}

    Color(double red, double green, double blue, double alpha) : Color() {
        r = red;
        g = green;
        b = blue;
        a = alpha;
    }

    static const Color White() { return Color(1.0, 1.0, 1.0); }

    static const Color Black() { return Color(0.0, 0.0, 0.0); }

    static const Color Gray() { return Color(0.5, 0.5, 0.5); }

    static const Color Red() { return Color(1.0, 0.0, 0.0); }

    static const Color Green() { return Color(0.0, 0.7, 0.0); }

    static const Color Blue() { return Color(0.0, 0.0, 1.0); }

    static const Color SteelBlue() { return Color(0.4, 0.7, 1.0); }

    static const Color Yellow() { return Color(1.0, 1.0, 0.0); }

    static const Color Orange() { return Color(1.0, 0.5, 0.0); }

    static const Color Purple() { return Color(0.5, 0.0, 1.0); }

    static const Color Chartreuse() { return Color(0.5, 1.0, 0.0); }

    static const Color Teal() { return Color(0.0, 1.0, 1.0); }

    static const Color Pink() { return Color(1.0, 0.0, 0.5); }
};


class VisualUtils {
private:
    static ros::Publisher nullpub;
public:

    static void PublishTrajToEllipsoid(ros::Publisher publisher_, Trajectory traj, int samples = 100){
        visualization_msgs::Marker ellipsoidMarker;
        visualization_msgs::MarkerArray ellipsoidMarkers;

        ellipsoidMarker.id = 0;
        ellipsoidMarker.type = visualization_msgs::Marker::MESH_RESOURCE;
        ellipsoidMarker.mesh_resource = "package://kino_planner/mesh/sphere.stl";
        ellipsoidMarker.header.stamp = ros::Time::now();
        ellipsoidMarker.header.frame_id = "world";
        ellipsoidMarker.pose.orientation.w = 1.00;
        ellipsoidMarker.action = visualization_msgs::Marker::ADD;
        ellipsoidMarker.ns = "ellipsoids";
        ellipsoidMarker.color.r = 1;
        ellipsoidMarker.color.g = 0;
        ellipsoidMarker.color.b = 1;
        ellipsoidMarker.color.a = 0.5;
        ellipsoidMarker.scale.x = 0.2 * 2.0;
        ellipsoidMarker.scale.y = 0.2 * 2.0;
        ellipsoidMarker.scale.z = 0.08 * 2.0;

        ellipsoidMarker.action = visualization_msgs::Marker::DELETEALL;
        ellipsoidMarkers.markers.push_back(ellipsoidMarker);
        publisher_.publish(ellipsoidMarkers);
        ellipsoidMarker.action = visualization_msgs::Marker::ADD;
        ellipsoidMarkers.markers.clear();

        double dt = traj.getTotalDuration() / samples;
        geometry_msgs::Point point;
        Eigen::Vector3d pos;
        Eigen::Matrix3d rotM;
        Eigen::Quaterniond quat;
        for (int i = 0; i <= samples; i++)
        {
            pos = traj.getPos(dt * i);
            ellipsoidMarker.pose.position.x = pos(0);
            ellipsoidMarker.pose.position.y = pos(1);
            ellipsoidMarker.pose.position.z = pos(2);
            traj.getRotation(dt * i, 0.0, 9.81, rotM);
            quat = Eigen::Quaterniond(rotM);
            ellipsoidMarker.pose.orientation.w = quat.w();
            ellipsoidMarker.pose.orientation.x = quat.x();
            ellipsoidMarker.pose.orientation.y = quat.y();
            ellipsoidMarker.pose.orientation.z = quat.z();
            ellipsoidMarkers.markers.push_back(ellipsoidMarker);
            ellipsoidMarker.id++;
        }

        publisher_.publish(ellipsoidMarkers);
    }

    static void PublishBallToMarkerArrays(ros::Publisher publisher_,Balls searched_pts){

        static int id = 0;
        visualization_msgs::MarkerArray mkr_arr;
        visualization_msgs::Marker  ball_list;
        ball_list.action = visualization_msgs::Marker::DELETEALL;
        mkr_arr.markers.push_back(ball_list);
        publisher_.publish(mkr_arr);
        mkr_arr.markers.clear();
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

            mkr_arr.markers.push_back(ball_list);
        }
        publisher_.publish(mkr_arr);
    }

    static void CalculatePoseFromAcc(Vec3 acc, geometry_msgs::Pose &out_pose) {
        Vec3 gravity(0, 0, 9.8);
        double yaw = 0;
        double a_T = (gravity + acc).norm();
        Eigen::Vector3d xB, yB, zB;
        Eigen::Vector3d xC(cos(yaw), sin(yaw), 0);

        zB = (gravity + acc).normalized();
        yB = ((zB).cross(xC)).normalized();
        xB = yB.cross(zB);
        Eigen::Matrix3d R;
        R << xB, yB, zB;
        Eigen::Quaterniond q(R);
        out_pose.orientation.w = q.w();
        out_pose.orientation.x = q.x();
        out_pose.orientation.y = q.y();
        out_pose.orientation.z = q.z();
    }

    static void PublishStartAndEndPointToMarkerArray(ros::Publisher &pub_, const geometry_msgs::Pose &start_,
                                                     const geometry_msgs::Pose &end_) {
        visualization_msgs::MarkerArray mkr_arr;
        double uav_scale = 0.003;
        Eigen::Vector3d uav_t_bias(0.5, -0.05, -0.32);
        Eigen::Quaterniond uav_q_bias(1, 0, 0, 0);
        {
            visualization_msgs::Marker marker;
            geometry_msgs::Pose pose = start_;
            Eigen::Quaterniond odom_q(pose.orientation.w, pose.orientation.x, pose.orientation.y, pose.orientation.z);
            Eigen::Vector3d odom_t(pose.position.x, pose.position.y, pose.position.z);
            Eigen::Matrix3d R = odom_q.toRotationMatrix() * uav_q_bias.toRotationMatrix();
            odom_q = Eigen::Quaterniond(R);
            odom_t +=   uav_t_bias ;


            // Set the frame ID and timestamp.  See the TF tutorials for information on these.
            marker.header.frame_id = "world";
            marker.header.stamp = ros::Time::now();

            // Set the namespace and id for this marker.  This serves to create a unique ID
            // Any marker sent with the same namespace and id will overwrite the old one
            marker.ns = "start";
            marker.id = 0;

            // Set the marker type
            marker.type = visualization_msgs::Marker::MESH_RESOURCE;

            // Set the marker action.  Options are ADD and DELETE
            marker.action = visualization_msgs::Marker::ADD;

            // Set the pose of the marker.  This is a full 6DOF pose relative to the frame/time specified in the header
            marker.pose.position.x = odom_t.x();
            marker.pose.position.y = odom_t.y();
            marker.pose.position.z = odom_t.z();
            marker.pose.orientation.w = odom_q.w();
            marker.pose.orientation.x = odom_q.x();
            marker.pose.orientation.y = odom_q.y();
            marker.pose.orientation.z = odom_q.z();


            // Set the scale of the marker -- 1x1x1 here means 1m on a side
            marker.scale.x = uav_scale;
            marker.scale.y = uav_scale;
            marker.scale.z = uav_scale;

            // Set the color -- be sure to set alpha to something non-zero!
            marker.color.r = 1;
            marker.color.g = 1;
            marker.color.b = 1;
            marker.color.a = 1.0;
            marker.lifetime = ros::Duration();

            // Publish the marker
            marker.mesh_resource = "package://kino_planner/mesh/hongyan.STL";
            mkr_arr.markers.push_back(marker);
        }

        {
            visualization_msgs::Marker marker;
            geometry_msgs::Pose pose = end_;
            Eigen::Quaterniond odom_q(pose.orientation.w, pose.orientation.x, pose.orientation.y, pose.orientation.z);
            Eigen::Vector3d odom_t(pose.position.x, pose.position.y, pose.position.z);
            Eigen::Matrix3d R = odom_q.toRotationMatrix() * uav_q_bias.toRotationMatrix();
            odom_q = Eigen::Quaterniond(R);
            odom_t += uav_t_bias;


            // Set the frame ID and timestamp.  See the TF tutorials for information on these.
            marker.header.frame_id = "world";
            marker.header.stamp = ros::Time::now();

            // Set the namespace and id for this marker.  This serves to create a unique ID
            // Any marker sent with the same namespace and id will overwrite the old one
            marker.ns = "end";
            marker.id = 1;

            // Set the marker type
            marker.type = visualization_msgs::Marker::MESH_RESOURCE;

            // Set the marker action.  Options are ADD and DELETE
            marker.action = visualization_msgs::Marker::ADD;

            // Set the pose of the marker.  This is a full 6DOF pose relative to the frame/time specified in the header
            marker.pose.position.x = odom_t.x();
            marker.pose.position.y = odom_t.y();
            marker.pose.position.z = odom_t.z();
            marker.pose.orientation.w = odom_q.w();
            marker.pose.orientation.x = odom_q.x();
            marker.pose.orientation.y = odom_q.y();
            marker.pose.orientation.z = odom_q.z();


            // Set the scale of the marker -- 1x1x1 here means 1m on a side
            marker.scale.x = uav_scale;
            marker.scale.y = uav_scale;
            marker.scale.z = uav_scale;

            // Set the color -- be sure to set alpha to something non-zero!
            marker.color.r = 1;
            marker.color.g = 1;
            marker.color.b = 1;
            marker.color.a = 1.0;
            marker.lifetime = ros::Duration();

            // Publish the marker
            marker.mesh_resource = "package://kino_planner/mesh/hongyan.STL";
            mkr_arr.markers.push_back(marker);
        }

        pub_.publish(mkr_arr);
    }

    static void VisualizeTopo(ros::Publisher pub_, const std::vector<Eigen::Vector3d> &p_head,
                              const std::vector<Eigen::Vector3d> &tracks,
                              ros::Time local_time) {
        if (tracks.empty() || p_head.empty() || tracks.size() != p_head.size()) {
            printf(" EMPTY TOPO\n");
            return;
        }
        visualization_msgs::MarkerArray mkrs;
        visualization_msgs::Marker topo;
        geometry_msgs::Point p;

        for (int i = 0; i < tracks.size(); i++) {
            p.x = p_head[i][0];
            p.y = p_head[i][1];
            p.z = p_head[i][2];
            topo.points.push_back(p);
            p.x += tracks[i][0];
            p.y += tracks[i][1];
            p.z += tracks[i][2];
            topo.points.push_back(p);
        }

        topo.header.frame_id = "world";
        topo.header.stamp = ros::Time::now();
        topo.ns = "topo";
        topo.action = visualization_msgs::Marker::ADD;
//        topo.lifetime = ros::Duration(0);
        topo.pose.orientation.w = 1.0;
        topo.id = 117;
        topo.type = visualization_msgs::Marker::LINE_LIST;
        topo.scale.x = 0.15;
        topo.color = Color::Green();
        topo.color.a = 1.0;
        mkrs.markers.push_back(topo);

        pub_.publish(mkrs);
    }

    static void PublishVectorVec3ToMarkerArray(ros::Publisher &publisher_, vector<Vec3> &points,
                                               Vec3 color = ORANGE, string pt_ns = "vecpts") {
        visualization_msgs::MarkerArray mrkarr;
        color /= 255;
        visualization_msgs::Marker mk;
//        mk.action = visualization_msgs::Marker::DELETEALL;
//        mrkarr.markers.push_back(mk);

        static int id = 0;
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

    static void
    PublishPieceToMarkerArray(ros::Publisher &pub_, Piece piece, string ns = "line", Vec3 color = Vec3(1, 0, 0)) {
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

    static void PublishFullTreeToMarkerArray(ros::Publisher &pub_, StateKdTree::Ptr skd_tree_) {
        double square_dis = 100 * 100;
        TreeNodePtr search_pt(new TreeNode);
        search_pt->x.col(0) = Vec3(0, 0, 0);
        vector<TreeNodePtr> out = skd_tree_->RadiusSearch(square_dis, search_pt);
        visualization_msgs::MarkerArray mrkarr;
        visualization_msgs::Marker mk;
        mk.action = visualization_msgs::Marker::DELETEALL;
        mrkarr.markers.push_back(mk);
        int cnt = 0;
        double ball_size = 0.2;
        int line_cnt = 0;
        print(fg(color::cornflower_blue), " -- [Visual] SKD size = {}\n", out.size());
        for (size_t i = 0; i < out.size(); i++) {

            TreeNodePtr node1 = out[i];
            TreeNodePtr node2 = out[i]->came_from;
            if (node2 == NULL) {
                continue;
            }
            visualization_msgs::Marker node_posi;
            node_posi.header.frame_id = "world";
            node_posi.header.stamp = ros::Time::now();
            node_posi.ns = "node";
            node_posi.id = cnt++;
            node_posi.action = visualization_msgs::Marker::ADD;
            node_posi.pose.orientation.w = 1.0;
            node_posi.type = visualization_msgs::Marker::SPHERE_LIST;
            // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
            node_posi.scale.x = ball_size;
            node_posi.scale.y = ball_size;
            node_posi.scale.z = ball_size;
            // Line list is blue
            node_posi.color.b = 0.5;
            node_posi.color.r = 0.5;
            node_posi.color.a = 1.0;
            // Create the vertices for the points and lines
            geometry_msgs::Point p;
            p.x = node1->x.col(0).x();
            p.y = node1->x.col(0).y();
            p.z = node1->x.col(0).z();
            node_posi.points.push_back(p);

            mrkarr.markers.push_back(node_posi);
            {
                PublishPieceToMarkerArray(pub_, node1->piece, "branch", SPRINGGREEN);
                visualization_msgs::Marker node_posi;
                node_posi.header.frame_id = "world";
                node_posi.header.stamp = ros::Time::now();
                node_posi.ns = "node";
                node_posi.id = cnt++;
                node_posi.action = visualization_msgs::Marker::ADD;
                node_posi.pose.orientation.w = 1.0;
                node_posi.type = visualization_msgs::Marker::SPHERE_LIST;
                // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                node_posi.scale.x = ball_size;
                node_posi.scale.y = ball_size;
                node_posi.scale.z = ball_size;
                // Line list is blue
                node_posi.color.b = 0.5;
                node_posi.color.r = 0.5;
                node_posi.color.a = 1.0;
                // Create the vertices for the points and lines
                geometry_msgs::Point p;
                p.x = node2->x.col(0).x();
                p.y = node2->x.col(0).y();
                p.z = node2->x.col(0).z();
                node_posi.points.push_back(p);
                mrkarr.markers.push_back(node_posi);

                // publish lines
                visualization_msgs::Marker line_list;
                line_list.header.frame_id = "world";
                line_list.header.stamp = ros::Time::now();
                line_list.ns = "lines";
                line_list.id = line_cnt++;
                line_list.action = visualization_msgs::Marker::ADD;
                line_list.pose.orientation.w = 1.0;
                line_list.type = visualization_msgs::Marker::LINE_LIST;
                // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                line_list.scale.x = 0.03;
                // Line list is blue
                line_list.color.b = 1.0;
                line_list.color.a = 1.0;
                // Create the vertices for the points and lines

                p.x = node1->x.col(0).x();
                p.y = node1->x.col(0).y();
                p.z = node1->x.col(0).z();
                // The line list needs two points for each line
                line_list.points.push_back(p);
                p.x = node2->x.col(0).x();
                p.y = node2->x.col(0).y();
                p.z = node2->x.col(0).z();
                line_list.points.push_back(p);

                mrkarr.markers.push_back(line_list);
            }
        }
        pub_.publish(mrkarr);
    }

    static void PublishNodeTreeToMarkerArray(ros::Publisher &pub_, StateKdTree::Ptr skd_tree_) {
        double square_dis = 100 * 100;
        TreeNodePtr search_pt(new TreeNode);
        search_pt->x.col(0) = Vec3(0, 0, 0);
        vector<TreeNodePtr> out = skd_tree_->RadiusSearch(square_dis, search_pt);
        visualization_msgs::MarkerArray mrkarr;
        visualization_msgs::Marker mk;
        mk.action = visualization_msgs::Marker::DELETEALL;
        mrkarr.markers.push_back(mk);
        int cnt = 0;
        double ball_size = 0.2;
        int line_cnt = 0;
        print(fg(color::cornflower_blue), " -- [Visual] SKD size = {}\n", out.size());
        for (size_t i = 0; i < out.size(); i++) {

            TreeNodePtr node1 = out[i];

            visualization_msgs::Marker node_posi;
            node_posi.header.frame_id = "world";
            node_posi.header.stamp = ros::Time::now();
            node_posi.ns = "node";
            node_posi.id = cnt++;
            node_posi.action = visualization_msgs::Marker::ADD;
            node_posi.pose.orientation.w = 1.0;
            node_posi.type = visualization_msgs::Marker::SPHERE_LIST;
            // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
            node_posi.scale.x = ball_size;
            node_posi.scale.y = ball_size;
            node_posi.scale.z = ball_size;
            // Line list is blue
            node_posi.color.b = 0.5;
            node_posi.color.r = 0.5;
            node_posi.color.a = 1.0;
            // Create the vertices for the points and lines
            geometry_msgs::Point p;
            p.x = node1->x.col(0).x();
            p.y = node1->x.col(0).y();
            p.z = node1->x.col(0).z();
            node_posi.points.push_back(p);

            mrkarr.markers.push_back(node_posi);

            TreeNodePtr node2 = out[i]->came_from;
            if (node2 != NULL) {

                visualization_msgs::Marker node_posi;
                node_posi.header.frame_id = "world";
                node_posi.header.stamp = ros::Time::now();
                node_posi.ns = "node";
                node_posi.id = cnt++;
                node_posi.action = visualization_msgs::Marker::ADD;
                node_posi.pose.orientation.w = 1.0;
                node_posi.type = visualization_msgs::Marker::SPHERE_LIST;
                // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                node_posi.scale.x = ball_size;
                node_posi.scale.y = ball_size;
                node_posi.scale.z = ball_size;
                // Line list is blue
                node_posi.color.b = 0.5;
                node_posi.color.r = 0.5;
                node_posi.color.a = 1.0;
                // Create the vertices for the points and lines
                geometry_msgs::Point p;
                p.x = node2->x.col(0).x();
                p.y = node2->x.col(0).y();
                p.z = node2->x.col(0).z();
                node_posi.points.push_back(p);
                mrkarr.markers.push_back(node_posi);

                // publish lines
                visualization_msgs::Marker line_list;
                line_list.header.frame_id = "world";
                line_list.header.stamp = ros::Time::now();
                line_list.ns = "lines";
                line_list.id = line_cnt++;
                line_list.action = visualization_msgs::Marker::ADD;
                line_list.pose.orientation.w = 1.0;
                line_list.type = visualization_msgs::Marker::LINE_LIST;
                // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                line_list.scale.x = 0.03;
                // Line list is blue
                line_list.color.b = 1.0;
                line_list.color.a = 1.0;
                // Create the vertices for the points and lines

                p.x = node1->x.col(0).x();
                p.y = node1->x.col(0).y();
                p.z = node1->x.col(0).z();
                // The line list needs two points for each line
                line_list.points.push_back(p);
                p.x = node2->x.col(0).x();
                p.y = node2->x.col(0).y();
                p.z = node2->x.col(0).z();
                line_list.points.push_back(p);

                mrkarr.markers.push_back(line_list);
            }
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

/*--------------  TRAJ  -------------------------------------*/
    static void publishTrajToPath(ros::Publisher &pub_, Trajectory traj) {
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


    static void PublishTrajPosiToMarkerArray(ros::Publisher &pub_, Trajectory traj, Vec3 color = YELLOW) {
        double t_sum = traj.getTotalDuration();
        double sum_dis = traj.getTotalLength();
        double eval_t = 0.01;
        double v_n;
        Vec3 cur_vel_dir, last_pos, cur_pos, end_point;
        Vec3 v2 = (traj.getPos(-1) - traj.getPos(0)).normalized();
        visualization_msgs::Marker line_list;
        line_list.action = visualization_msgs::Marker::DELETEALL;
        visualization_msgs::MarkerArray mrkarr;
        mrkarr.markers.push_back(line_list);
        pub_.publish(mrkarr);
        color /= 255;
        mrkarr.markers.clear();
        int idx = 0;
        cur_pos = traj.getPos(0);
        while (eval_t + 1e-4 < t_sum) {
            visualization_msgs::Marker line_list;
            last_pos = cur_pos;
            cur_pos = traj.getPos(eval_t);

            {
                line_list.header.frame_id = "world";
                line_list.header.stamp = ros::Time::now();
                line_list.ns = "posi";
                line_list.id = idx++;
                line_list.action = visualization_msgs::Marker::ADD;
                line_list.pose.orientation.w = 1.0;
                line_list.type = visualization_msgs::Marker::LINE_LIST;
                // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
                line_list.scale.x = 0.1;
                // Line list is red
                line_list.color.r = color[0];
                line_list.color.g = color[1];
                line_list.color.b = color[2];
                line_list.color.a = 1;
                // Create the vertices for the points and lines
                geometry_msgs::Point p;
                p.x = last_pos.x();
                p.y = last_pos.y();
                p.z = last_pos.z();
                // The line list needs two points for each line
                line_list.points.push_back(p);
                p.x = cur_pos.x();
                p.y = cur_pos.y();
                p.z = cur_pos.z();
                line_list.points.push_back(p);

                mrkarr.markers.push_back(line_list);
            }
            eval_t += t_sum / sum_dis / 32;

        }
        pub_.publish(mrkarr);
    }

    static void PublishTrajVelToMarkerArray(ros::Publisher &pub_, Trajectory traj, Vec3 color = RED) {
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

    static void PublishTrajAccToMarkerArray(ros::Publisher &pub_, Trajectory traj, Vec3 color = BLUE) {
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

    static void PublishTrajJerkToMarkerArray(ros::Publisher &pub_, Trajectory traj, Vec3 color = BLUE) {
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
            cur_acc = traj.getJerk(eval_t);
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

    static void PublishTrajSnapToMarkerArray(ros::Publisher &pub_, Trajectory traj, Vec3 color = GREENYELLOW) {
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
    }

};


#endif //SRC_KRRT_PLANNER_H
