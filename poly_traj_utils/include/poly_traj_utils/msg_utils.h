//
// Created by yunfan on 2021/8/18.
//

#ifndef SRC_MSG_UTILS_H
#define SRC_MSG_UTILS_H
#include "Eigen/Dense"
#include "tf/tf.h"
#include "fmt/color.h"

using namespace std;
using namespace Eigen;
typedef Matrix<double, 3, 3> Mat33;
typedef Matrix<double, 3, 1> Vec3;
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
#endif //SRC_MSG_UTILS_H
