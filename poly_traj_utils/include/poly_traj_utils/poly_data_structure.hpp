/*
    MIT License

    Copyright (c) 2021 Zhepei Wang (wangzhepei@live.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

/*
    @article{WANG2021GCOPTER,
        title={Geometrically Constrained Trajectory Optimization for Multicopters},
        author={Wang, Zhepei and Zhou, Xin and Xu, Chao and Gao, Fei},
        journal={arXiv preprint arXiv:2103.00190},
        year={2021}
    }
*/

#ifndef POLY_DATA_STRUCTURE
#define POLY_DATA_STRUCTURE

#include "root_finder.hpp"
#include <vector>
#include <list>
#include <Eigen/Eigen>
#include <iostream>
#include "fmt/color.h"
#include "scope_timer.hpp"

using namespace std;


// Type for piece boundary condition and coefficient matrix
//typedef Eigen::Matrix<double, TrajDim, TrajOrder + 1> BoundaryCond;
//typedef Eigen::Matrix<double, TrajDim, TrajOrder + 1> CoefficientMat;
//typedef Eigen::Matrix<double, TrajDim, TrajOrder> VelCoefficientMat;
//typedef Eigen::Matrix<double, TrajDim, TrajOrder - 1> AccCoefficientMat;
typedef Eigen::Matrix<double, 3, 1> Vec3;
typedef Eigen::Matrix<double, 3, 2> StatePV;
typedef Eigen::Matrix<double, 3, 3> StatePVA;
typedef Eigen::Matrix<double, 3, 4> StatePVAJ;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DynamicMat;

// A single piece of a trajectory, which is indeed a polynomial
class Piece {
private:
    // Piece(t) = c5*t^5 + c4*t^4 + ... + c1*t + c0
    // The natural coefficient matrix = [c5,c4,c3,c2,c1,c0]
    double duration;
    // Any time in [0, T] is normalized into [0.0, 1.0]
    // Therefore, nCoeffMat = [c5*T^5,c4*T^4,c3*T^3,c2*T^2,c1*T,c0*1]
    // is used for better numerical stability
    //    CoefficientMat nCoeffMat;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nCoeffMat;
    int order_;
    bool had_sampled_ = false;
    std::vector<Vec3> sampled_positions_;

    bool had_length_ = false;
    double length_;
    bool is_empty_{true};

    double cost_;

    DynamicMat fw_mat_;
    bool has_fw_mat_ = false;

public:
    Piece() = default;

    // Constructor from duration and coefficient
    Piece(double dur, DynamicMat coeffs) : duration(dur) {
        order_ = coeffs.cols() - 1;
        nCoeffMat.resize(3, coeffs.cols());
        double t = 1.0;
        for (int i = order_; i >= 0; i--) {
            nCoeffMat.col(i) = coeffs.col(i) * t;
            t *= dur;
        }
        is_empty_ = false;
    }

    // Constructor from boundary condition and duration
    Piece(DynamicMat boundCond, double dur) : duration(dur) {

        if (boundCond.cols() == 8) {
            order_ = 7;
            nCoeffMat.resize(3, order_ + 1);
            // The BoundaryCond matrix boundCond = [p(0),v(0),a(0),p(T),v(T),a(T)]
            /*                                          0   1    2    3    4   5     6    7
             *   The BoundaryCond matrix boundCond = [p(0),v(0),a(0),j(0),p(T),v(T),a(T),j(T)]
             * */
            double t1 = dur;
            double t2 = t1 * t1;
            double t3 = t2 * t1;


            // Inverse mapping is computed without explicit matrix inverse
            // It maps boundary condition to normalized coefficient matrix
            Eigen::Array3d p_0 = boundCond.col(0),
                    v_0 = boundCond.col(1),
                    a_0 = boundCond.col(2),
                    j_0 = boundCond.col(3),
                    p_f = boundCond.col(4),
                    v_f = boundCond.col(5),
                    a_f = boundCond.col(6),
                    j_f = boundCond.col(7);

            nCoeffMat.col(0) =
                    20 * p_0 - 20 * p_f + 10 * t1 * v_0 + 10 * t1 * v_f + 2 * t2 * a_0 - 2 * t2 * a_f + (t3 * j_0) / 6 +
                    (t3 * j_f) / 6;
            nCoeffMat.col(1) =
                    70 * p_f - 70 * p_0 - 36 * t1 * v_0 - 34 * t1 * v_f - (15 * t2 * a_0) / 2 + (13 * t2 * a_f) / 2 -
                    (2 * t3 * j_0) / 3 - (t3 * j_f) / 2;
            nCoeffMat.col(2) =
                    84 * p_0 - 84 * p_f + 45 * t1 * v_0 + 39 * t1 * v_f + 10 * t2 * a_0 - 7 * t2 * a_f + t3 * j_0 +
                    (t3 * j_f) / 2;
            nCoeffMat.col(3) = 35 * p_f - 35 * p_0 - 20 * t1 * v_0 - 15 * t1 * v_f - 5 * t2 * a_0 + (5 * t2 * a_f) / 2 -
                               (2 * t3 * j_0) / 3 - (t3 * j_f) / 6;
            nCoeffMat.col(4) = (t3 * j_0) / 6;
            nCoeffMat.col(5) = (t2 * a_0) / 2;
            nCoeffMat.col(6) = t1 * v_0;
            nCoeffMat.col(7) = p_0;
        } else if (boundCond.cols() == 6) {
            // The BoundaryCond matrix boundCond = [p(0),v(0),a(0),p(T),v(T),a(T)]
            order_ = 5;
            nCoeffMat.resize(3, order_ + 1);
            double t1 = dur;
            double t2 = t1 * t1;

            // Inverse mapping is computed without explicit matrix inverse
            // It maps boundary condition to normalized coefficient matrix
            nCoeffMat.col(0) = 0.5 * (boundCond.col(5) - boundCond.col(2)) * t2 -
                               3.0 * (boundCond.col(1) + boundCond.col(4)) * t1 +
                               6.0 * (boundCond.col(3) - boundCond.col(0));
            nCoeffMat.col(1) = (-boundCond.col(5) + 1.5 * boundCond.col(2)) * t2 +
                               (8.0 * boundCond.col(1) + 7.0 * boundCond.col(4)) * t1 +
                               15.0 * (-boundCond.col(3) + boundCond.col(0));
            nCoeffMat.col(2) = (0.5 * boundCond.col(5) - 1.5 * boundCond.col(2)) * t2 -
                               (6.0 * boundCond.col(1) + 4.0 * boundCond.col(4)) * t1 +
                               10.0 * (boundCond.col(3) - boundCond.col(0));
            nCoeffMat.col(3) = 0.5 * boundCond.col(2) * t2;
            nCoeffMat.col(4) = boundCond.col(1) * t1;
            nCoeffMat.col(5) = boundCond.col(0);
        }
        is_empty_ = false;
    }

    inline bool empty() {
        return is_empty_;
    }

    inline void setCost(double cost) {
        cost_ = cost;
    }

    inline double getCost() {
        return cost_;
    };

    inline void reset() {
        nCoeffMat.setZero();
        duration = 0;
        had_length_ = false;
        had_sampled_ = false;
        sampled_positions_.clear();
    }

    inline void resetDuration(double dur) {
        DynamicMat mat = getCoeffMat();
        duration = dur;
        double t = 1.0;
        cost_ = -1;
        had_sampled_ = false;
        for (int i = order_; i >= 0; i--) {
            nCoeffMat.col(i) = mat.col(i) * t;
            t *= dur;
        }
    }

    inline int getOrder() const {
        return order_;
    }

    inline double getDuration() const {
        return duration;
    }

    // Get the position at time t in this piece
    inline Eigen::Vector3d getPos(double t) const {
        if (t == -1) {
            t = duration;
        }
        // Normalize the time
        t /= duration;
        Eigen::Vector3d pos(0.0, 0.0, 0.0);
        double tn = 1.0;
        for (int i = order_; i >= 0; i--) {
            pos += tn * nCoeffMat.col(i);
            tn *= t;
        }
        // The pos is not affected by normalization
        return pos;
    }

    // Get the velocity at time t in this piece
    inline Eigen::Vector3d getVel(double t) const {
        if (t == -1) {
            t = duration;
        }
        // Normalize the time
        t /= duration;
        Eigen::Vector3d vel(0.0, 0.0, 0.0);
        double tn = 1.0;
        int n = 1;
        for (int i = order_ - 1; i >= 0; i--) {
            vel += n * tn * nCoeffMat.col(i);
            tn *= t;
            n++;
        }
        // Recover the actual vel
        vel /= duration;
        return vel;
    }

    // Get the acceleration at time t in this piece
    inline Eigen::Vector3d getAcc(double t) const {
        if (t == -1) {
            t = duration;
        }
        // Normalize the time
        t /= duration;
        Eigen::Vector3d acc(0.0, 0.0, 0.0);
        double tn = 1.0;
        int m = 1;
        int n = 2;
        for (int i = order_ - 2; i >= 0; i--) {
            acc += m * n * tn * nCoeffMat.col(i);
            tn *= t;
            m++;
            n++;
        }
        // Recover the actual acc
        acc /= duration * duration;
        return acc;
    }

    // Get the jerk at time t in this piece
    inline Eigen::Vector3d getJerk(double t) const {
        if (t == -1) {
            t = duration;
        }
        // Normalize the time
        t /= duration;
        Eigen::Vector3d jerk(0.0, 0.0, 0.0);
        double tn = 1.0;
        int m = 1;
        int n = 2;
        int k = 3;
        for (int i = order_ - 3; i >= 0; i--) {
            jerk += k * m * n * tn * nCoeffMat.col(i);
            tn *= t;
            k++;
            m++;
            n++;
        }
        // Recover the actual acc
        jerk /= duration * duration * duration;
        return jerk;
    }

    // Get the snap at time t in this piece
    inline Eigen::Vector3d getSnap(double t) const {
        if (t == -1) {
            t = duration;
        }
        // Normalize the time
        t /= duration;
        Eigen::Vector3d snap(0.0, 0.0, 0.0);
        double tn = 1.0;
        int m = 1;
        int n = 2;
        int k = 3;
        int w = 4;
        for (int i = order_ - 4; i >= 0; i--) {
            snap += w * k * m * n * tn * nCoeffMat.col(i);
            tn *= t;
            w++;
            k++;
            m++;
            n++;
        }
        // Recover the actual acc
        snap /= duration * duration * duration * duration;
        return snap;
    }

    // Get the boundary condition of this piece
    inline DynamicMat getBoundCond() const {
        DynamicMat boundCond;
        if (order_ == 7) {
            boundCond.resize(3, order_ + 1);
            boundCond << getPos(0.0), getVel(0.0), getAcc(0.0), getJerk(0.0),
                    getPos(duration), getVel(duration), getAcc(duration), getJerk(duration);
        } else if (order_ == 5) {
            boundCond.resize(3, order_ + 1);
            boundCond << getPos(0.0), getVel(0.0), getAcc(0.0),
                    getPos(duration), getVel(duration), getAcc(duration);
        }
        return boundCond;
    }

    // Get the coefficient matrix of the piece
    // Default arg chooses the natural coefficients
    // If normalized version is needed, set the arg true
    inline DynamicMat getCoeffMat(bool normalized = false) const {
        DynamicMat posCoeffsMat;
        posCoeffsMat.resize(3, order_ + 1);
        double t = 1;
        for (int i = order_; i >= 0; i--) {
            posCoeffsMat.col(i) = nCoeffMat.col(i) / t;
            t *= normalized ? 1.0 : duration;
        }
        return posCoeffsMat;
    }

    // Get the polynomial coefficients of velocity of this piece
    // Default arg chooses the natural coefficients
    // If normalized version is needed, set the arg true
    inline DynamicMat getVelCoeffMat(bool normalized = false) const {
        DynamicMat velCoeffMat;
        velCoeffMat.resize(3, order_);
        int n = 1;
        double t = 1.0;
        t *= normalized ? 1.0 : duration;
        for (int i = order_ - 1; i >= 0; i--) {
            velCoeffMat.col(i) = n * nCoeffMat.col(i) / t;
            n++;
            t *= normalized ? 1.0 : duration;
        }
        return velCoeffMat;
    }

    // Get the polynomial coefficients of acceleration of this piece
    // Default arg chooses the natural coefficients
    // If normalized version is needed, set the arg true
    inline DynamicMat getAccCoeffMat(bool normalized = false) const {
        DynamicMat accCoeffMat;
        accCoeffMat.resize(3, order_ - 1);

        int n = 2;
        int m = 1;
        double t = 1.0;
        t *= normalized ? 1.0 : duration * duration;
        for (int i = order_ - 2; i >= 0; i--) {
            accCoeffMat.col(i) = n * m * nCoeffMat.col(i) / t;
            n++;
            m++;
            t *= normalized ? 1.0 : duration;
        }
        return accCoeffMat;
    }

    // Get the max velocity rate of the piece
    inline double getMaxVelRate() const {
        // Compute normalized squared vel norm polynomial coefficient matrix
        Eigen::MatrixXd nVelCoeffMat = getVelCoeffMat(true);
        Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                                RootFinder::polySqr(nVelCoeffMat.row(1)) +
                                RootFinder::polySqr(nVelCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++) {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON) {
            return 0.0;
        } else {
            // Search an open interval whose boundaries are not zeros
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON) {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON) {
                r = 0.5 * (r + 1.0);
            }
            // Find all stationaries
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);

            // Check boundary points and stationaries within duration
            candidates.insert(0.0);
            candidates.insert(1.0);
            double maxVelRateSqr = -INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
                 it != candidates.end();
                 it++) {
                if (0.0 <= *it && 1.0 >= *it) {
                    // Recover the actual time then get the vel squared norm
                    tempNormSqr = getVel((*it) * duration).squaredNorm();
                    maxVelRateSqr = maxVelRateSqr < tempNormSqr ? tempNormSqr : maxVelRateSqr;
                }
            }
            return sqrt(maxVelRateSqr);
        }
    }

    // Get the max velocity rate of the piece
    inline double getMinVelRate() const {
        // Compute normalized squared vel norm polynomial coefficient matrix
        Eigen::MatrixXd nVelCoeffMat = getVelCoeffMat(true);
        Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                                RootFinder::polySqr(nVelCoeffMat.row(1)) +
                                RootFinder::polySqr(nVelCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++) {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON) {
            return 0.0;
        } else {
            // Search an open interval whose boundaries are not zeros
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON) {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON) {
                r = 0.5 * (r + 1.0);
            }
            // Find all stationaries
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);

            // Check boundary points and stationaries within duration
            candidates.insert(0.0);
            candidates.insert(1.0);
            double minVelRateSqr = INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
                 it != candidates.end();
                 it++) {
                if (0.0 <= *it && 1.0 >= *it) {
                    double cur_t = (*it) * duration;
                    if (cur_t < 0.01 || cur_t > duration - 0.01)
                        continue;
                    // Recover the actual time then get the vel squared norm
                    tempNormSqr = getVel(cur_t).squaredNorm();
                    minVelRateSqr = minVelRateSqr > tempNormSqr ? tempNormSqr : minVelRateSqr;
                }
            }
            return sqrt(minVelRateSqr);
        }
    }

    // Get the max acceleration rate of the piece
    inline double getMaxAccRate() const {
        // Compute normalized squared acc norm polynomial coefficient matrix
        Eigen::MatrixXd nAccCoeffMat = getAccCoeffMat(true);
        Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                                RootFinder::polySqr(nAccCoeffMat.row(1)) +
                                RootFinder::polySqr(nAccCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++) {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON) {
            return 0.0;
        } else {
            // Search an open interval whose boundaries are not zeros
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON) {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON) {
                r = 0.5 * (r + 1.0);
            }
            // Find all stationaries
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);
            // Check boundary points and stationaries within duration
            candidates.insert(0.0);
            candidates.insert(1.0);
            double maxAccRateSqr = -INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
                 it != candidates.end();
                 it++) {
                if (0.0 <= *it && 1.0 >= *it) {
                    // Recover the actual time then get the acc squared norm
                    tempNormSqr = getAcc((*it) * duration).squaredNorm();
                    maxAccRateSqr = maxAccRateSqr < tempNormSqr ? tempNormSqr : maxAccRateSqr;
                }
            }
            return sqrt(maxAccRateSqr);
        }
    }

    // Check whether velocity rate of the piece is always less than maxVelRate
    inline bool checkMaxVelRate(double maxVelRate) const {
        double sqrMaxVelRate = maxVelRate * maxVelRate;
        if (getVel(0.0).squaredNorm() >= sqrMaxVelRate ||
            getVel(duration).squaredNorm() >= sqrMaxVelRate) {
            return false;
        } else {
            Eigen::MatrixXd nVelCoeffMat = getVelCoeffMat(true);
            Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                                    RootFinder::polySqr(nVelCoeffMat.row(1)) +
                                    RootFinder::polySqr(nVelCoeffMat.row(2));
            // Convert the actual squared maxVelRate to a normalized one
            double t2 = duration * duration;
            coeff.tail<1>()(0) -= sqrMaxVelRate * t2;
            // Directly check the root existence in the normalized interval
            return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
        }
    }

    // Check whether accleration rate of the piece is always less than maxAccRate
    inline bool checkMaxAccRate(double maxAccRate) const {
        double sqrMaxAccRate = maxAccRate * maxAccRate;
        if (getAcc(0.0).squaredNorm() >= sqrMaxAccRate ||
            getAcc(duration).squaredNorm() >= sqrMaxAccRate) {
            return false;
        } else {
            Eigen::MatrixXd nAccCoeffMat = getAccCoeffMat(true);
            Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                                    RootFinder::polySqr(nAccCoeffMat.row(1)) +
                                    RootFinder::polySqr(nAccCoeffMat.row(2));
            // Convert the actual squared maxAccRate to a normalized one
            double t2 = duration * duration;
            double t4 = t2 * t2;
            coeff.tail<1>()(0) -= sqrMaxAccRate * t4;
            // Directly check the root existence in the normalized interval
            return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
        }
    }

    //Scale the Piece(t) to Piece(k*t)
    inline void scaleTime(double k) {
        duration /= k;
        return;
    }

    inline std::vector<Vec3> getTraj(double dt) {
        if (had_sampled_) {
            return sampled_positions_;
        }
        sampled_positions_.clear();
        had_sampled_ = true;
        for (double t = 0.0; t < duration; t += dt) {
            Eigen::Vector3d pos;
            sampled_positions_.push_back(getPos(t));
        }

        return sampled_positions_;
    }

    inline double getLength() {
        if (had_length_) {
            return length_;
        }
        length_ = 0;
        had_length_ = true;
        Vec3 pos = getPos(0), last_pos;
        for (double t = 0.0; t < duration; t += 0.01) {
            last_pos = pos;
            pos = getPos(t);
            length_ += (pos - last_pos).norm();
        }
        return length_;
    }


    /// Return: false means search processes occur error
    enum FW_RETCODE {
        ERROR = -1,
        SUCCESS = 1,
        END_IN_DIS = 2,
    };

    /// Do forward search in [search_time, duration] all t out of this range is blocked.
    inline int
    getForwardPosition(const double search_time, ///<[in] input search time on the piece, should be [0, duration]'
                       const double dist,    ///<[in
                       double &out_time) {

        /*   /// Typical usage:
         *
             double fw_t;
             int ret_code =  piece.getForwardPosition(search_time, cur_dist, fw_t);
             if(ret_code < 0)
             {
                 print(fg(color::red), " -- [CHK] Fowrdsearch error.\n");
             }
             next_pt = piece.getPos(fw_t);
             search_time = fw_t;
        */

        if(search_time > duration){
            fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, " -- [FWS] Forward search time out of traj.\n");
            return ERROR;
        }

        Vec3 bias_pt = getPos(search_time);
        Eigen::VectorXd coeffsGradT;
        std::set<double> roots;

        if (order_ == 5) {

            coeffsGradT.resize(11);
            DynamicMat coeff_mat = getCoeffMat(), c;
            c.resize(3, 6);

            for (int i = 0; i < 6; i++) {
                c.row(0)[i] = coeff_mat.row(0)[5 - i];
                c.row(1)[i] = coeff_mat.row(1)[5 - i];
                c.row(2)[i] = coeff_mat.row(2)[5 - i];
            }
            c.row(0)[0] -= bias_pt[0];
            c.row(1)[0] -= bias_pt[1];
            c.row(2)[0] -= bias_pt[2];

            coeffsGradT(0) = c(0, 5) * c(0, 5) + c(1, 5) * c(1, 5) + c(2, 5) * c(2, 5);
            coeffsGradT(1) = 2 * c(0, 4) * c(0, 5) + 2 * c(1, 4) * c(1, 5) + 2 * c(2, 4) * c(2, 5);
            coeffsGradT(2) = c(0, 4) * c(0, 4) + c(1, 4) * c(1, 4) + c(2, 4) * c(2, 4) + 2 * c(0, 3) * c(0, 5) +
                             2 * c(1, 3) * c(1, 5) + 2 * c(2, 3) * c(2, 5);
            coeffsGradT(3) =
                    2 * c(0, 2) * c(0, 5) + 2 * c(0, 3) * c(0, 4) + 2 * c(1, 2) * c(1, 5) + 2 * c(1, 3) * c(1, 4) +
                    2 * c(2, 2) * c(2, 5) + 2 * c(2, 3) * c(2, 4);
            coeffsGradT(4) = c(0, 3) * c(0, 3) + c(1, 3) * c(1, 3) + c(2, 3) * c(2, 3) + 2 * c(0, 1) * c(0, 5) +
                             2 * c(0, 2) * c(0, 4) + 2 * c(1, 1) * c(1, 5) + 2 * c(1, 2) * c(1, 4) +
                             2 * c(2, 1) * c(2, 5) + 2 * c(2, 2) * c(2, 4);
            coeffsGradT(5) =
                    2 * c(0, 0) * c(0, 5) + 2 * c(0, 1) * c(0, 4) + 2 * c(0, 2) * c(0, 3) + 2 * c(1, 0) * c(1, 5) +
                    2 * c(1, 1) * c(1, 4) + 2 * c(1, 2) * c(1, 3) + 2 * c(2, 0) * c(2, 5) + 2 * c(2, 1) * c(2, 4) +
                    2 * c(2, 2) * c(2, 3);
            coeffsGradT(6) = c(0, 2) * c(0, 2) + c(1, 2) * c(1, 2) + c(2, 2) * c(2, 2) + 2 * c(0, 0) * c(0, 4) +
                             2 * c(0, 1) * c(0, 3) + 2 * c(1, 0) * c(1, 4) + 2 * c(1, 1) * c(1, 3) +
                             2 * c(2, 0) * c(2, 4) + 2 * c(2, 1) * c(2, 3);
            coeffsGradT(7) =
                    2 * c(0, 0) * c(0, 3) + 2 * c(0, 1) * c(0, 2) + 2 * c(1, 0) * c(1, 3) + 2 * c(1, 1) * c(1, 2) +
                    2 * c(2, 0) * c(2, 3) + 2 * c(2, 1) * c(2, 2);
            coeffsGradT(8) = c(0, 1) * c(0, 1) + c(1, 1) * c(1, 1) + c(2, 1) * c(2, 1) + 2 * c(0, 0) * c(0, 2) +
                             2 * c(1, 0) * c(1, 2) + 2 * c(2, 0) * c(2, 2);
            coeffsGradT(9) = 2 * c(0, 0) * c(0, 1) + 2 * c(1, 0) * c(1, 1) + 2 * c(2, 0) * c(2, 1);
            coeffsGradT(10) = c(0, 0) * c(0, 0) + c(1, 0) * c(1, 0) + c(2, 0) * c(2, 0);
            coeffsGradT(10) -= dist * dist;

            roots = RootFinder::solvePolynomial(coeffsGradT, search_time, duration, 1e-3);

        } else if (order_ == 7) {

            coeffsGradT.resize(15);
            DynamicMat coeff_mat = getCoeffMat(), c;
            c.resize(3, 8);
            Vec3 bias_pt = getPos(search_time);
            for (int i = 0; i < 8; i++) {
                c.row(0)[i] = coeff_mat.row(0)[7 - i];
                c.row(1)[i] = coeff_mat.row(1)[7 - i];
                c.row(2)[i] = coeff_mat.row(2)[7 - i];
            }
            c.row(0)[0] -= bias_pt[0];
            c.row(1)[0] -= bias_pt[1];
            c.row(2)[0] -= bias_pt[2];

            coeffsGradT(0) = c(0, 7) * c(0, 7) + c(1, 7) * c(1, 7) + c(2, 7) * c(2, 7);
            coeffsGradT(1) = 2 * c(0, 6) * c(0, 7) + 2 * c(1, 6) * c(1, 7) + 2 * c(2, 6) * c(2, 7);
            coeffsGradT(2) = c(0, 6) * c(0, 6) + c(1, 6) * c(1, 6) + c(2, 6) * c(2, 6) + 2 * c(0, 5) * c(0, 7) +
                             2 * c(1, 5) * c(1, 7) + 2 * c(2, 5) * c(2, 7);
            coeffsGradT(3) =
                    2 * c(0, 4) * c(0, 7) + 2 * c(0, 5) * c(0, 6) + 2 * c(1, 4) * c(1, 7) + 2 * c(1, 5) * c(1, 6) +
                    2 * c(2, 4) * c(2, 7) + 2 * c(2, 5) * c(2, 6);
            coeffsGradT(4) = c(0, 5) * c(0, 5) + c(1, 5) * c(1, 5) + c(2, 5) * c(2, 5) + 2 * c(0, 3) * c(0, 7) +
                             2 * c(0, 4) * c(0, 6) + 2 * c(1, 3) * c(1, 7) + 2 * c(1, 4) * c(1, 6) +
                             2 * c(2, 3) * c(2, 7) + 2 * c(2, 4) * c(2, 6);
            coeffsGradT(5) =
                    2 * c(0, 2) * c(0, 7) + 2 * c(0, 3) * c(0, 6) + 2 * c(0, 4) * c(0, 5) + 2 * c(1, 2) * c(1, 7) +
                    2 * c(1, 3) * c(1, 6) + 2 * c(1, 4) * c(1, 5) + 2 * c(2, 2) * c(2, 7) + 2 * c(2, 3) * c(2, 6) +
                    2 * c(2, 4) * c(2, 5);
            coeffsGradT(6) = c(0, 4) * c(0, 4) + c(1, 4) * c(1, 4) + c(2, 4) * c(2, 4) + 2 * c(0, 1) * c(0, 7) +
                             2 * c(0, 2) * c(0, 6) + 2 * c(0, 3) * c(0, 5) + 2 * c(1, 1) * c(1, 7) +
                             2 * c(1, 2) * c(1, 6) + 2 * c(1, 3) * c(1, 5) + 2 * c(2, 1) * c(2, 7) +
                             2 * c(2, 2) * c(2, 6) + 2 * c(2, 3) * c(2, 5);
            coeffsGradT(7) =
                    2 * c(0, 0) * c(0, 7) + 2 * c(0, 1) * c(0, 6) + 2 * c(0, 2) * c(0, 5) + 2 * c(0, 3) * c(0, 4) +
                    2 * c(1, 0) * c(1, 7) + 2 * c(1, 1) * c(1, 6) + 2 * c(1, 2) * c(1, 5) + 2 * c(1, 3) * c(1, 4) +
                    2 * c(2, 0) * c(2, 7) + 2 * c(2, 1) * c(2, 6) + 2 * c(2, 2) * c(2, 5) + 2 * c(2, 3) * c(2, 4);
            coeffsGradT(8) = c(0, 3) * c(0, 3) + c(1, 3) * c(1, 3) + c(2, 3) * c(2, 3) + 2 * c(0, 0) * c(0, 6) +
                             2 * c(0, 1) * c(0, 5) + 2 * c(0, 2) * c(0, 4) + 2 * c(1, 0) * c(1, 6) +
                             2 * c(1, 1) * c(1, 5) + 2 * c(1, 2) * c(1, 4) + 2 * c(2, 0) * c(2, 6) +
                             2 * c(2, 1) * c(2, 5) + 2 * c(2, 2) * c(2, 4);
            coeffsGradT(9) =
                    2 * c(0, 0) * c(0, 5) + 2 * c(0, 1) * c(0, 4) + 2 * c(0, 2) * c(0, 3) + 2 * c(1, 0) * c(1, 5) +
                    2 * c(1, 1) * c(1, 4) + 2 * c(1, 2) * c(1, 3) + 2 * c(2, 0) * c(2, 5) + 2 * c(2, 1) * c(2, 4) +
                    2 * c(2, 2) * c(2, 3);
            coeffsGradT(10) = c(0, 2) * c(0, 2) + c(1, 2) * c(1, 2) + c(2, 2) * c(2, 2) + 2 * c(0, 0) * c(0, 4) +
                              2 * c(0, 1) * c(0, 3) + 2 * c(1, 0) * c(1, 4) + 2 * c(1, 1) * c(1, 3) +
                              2 * c(2, 0) * c(2, 4) + 2 * c(2, 1) * c(2, 3);
            coeffsGradT(11) =
                    2 * c(0, 0) * c(0, 3) + 2 * c(0, 1) * c(0, 2) + 2 * c(1, 0) * c(1, 3) + 2 * c(1, 1) * c(1, 2) +
                    2 * c(2, 0) * c(2, 3) + 2 * c(2, 1) * c(2, 2);
            coeffsGradT(12) = c(0, 1) * c(0, 1) + c(1, 1) * c(1, 1) + c(2, 1) * c(2, 1) + 2 * c(0, 0) * c(0, 2) +
                              2 * c(1, 0) * c(1, 2) + 2 * c(2, 0) * c(2, 2);
            coeffsGradT(13) = 2 * c(0, 0) * c(0, 1) + 2 * c(1, 0) * c(1, 1) + 2 * c(2, 0) * c(2, 1);
            coeffsGradT(14) = c(0, 0) * c(0, 0) + c(1, 0) * c(1, 0) + c(2, 0) * c(2, 0);
            coeffsGradT(14) -= dist * dist;
            roots = RootFinder::solvePolynomial(coeffsGradT, search_time, duration, 1e-3);
        } else {
            fmt::print(fg(fmt::color::red), " -- [PIECE] Wrong piece order, force return.\n");
            return ERROR;
        }


        if (roots.size() == 0) {
            double end_dis = (bias_pt - getPos(-1)).norm();
            if (end_dis < dist + 0.01) {
                out_time = duration;
                return END_IN_DIS;
            }
            fmt::print(" -- [FWS] Cur t = {}, total_t = {} , forward_dis = {}, end_dis = {}\n", search_time, duration, dist,
                       end_dis);
            fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, " -- [FWS] Cannot get forward search.\n");
            return ERROR;
        }

        for (const double &root : roots) {
            if (root > duration) {
                continue;
            }
            if (root - search_time >= 0) {
                out_time = root;
                return SUCCESS;
            }
        }
        fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "All roots larger than duration({}).\n", duration);
        return ERROR;
    }

};

// A whole trajectory which contains multiple pieces
class Trajectory {
private:
    typedef std::vector<Piece> Pieces;
    Pieces pieces;
    bool had_sampled_ = false;
    std::vector<Vec3> sampled_positions_;

    bool had_length_ = false;
    double length_;
    double cost_;
public:
    Trajectory() = default;

    // Constructor from durations and coefficient matrices
    Trajectory(const std::vector<double> &durs,
               const std::vector<DynamicMat> &coeffMats) {
        int N = std::min(durs.size(), coeffMats.size());
        for (int i = 0; i < N; i++) {
            pieces.emplace_back(durs[i], coeffMats[i]);
        }
    }

    inline void setCost(double cost) {
        cost_ = cost;
    }

    inline double getCost() {
        return cost_;
    }

    inline size_t getPieceNum() const {
        return pieces.size();
    }

    inline void reset() {
        pieces.clear();
        had_length_ = false;
        had_sampled_ = false;
        sampled_positions_.clear();
    }

    // Get durations vector of all pieces
    inline std::vector<double> getDurations() const {
        std::vector<double> durations;
        durations.reserve(getPieceNum());
        for (int i = 0; i < getPieceNum(); i++) {
            durations.push_back(pieces[i].getDuration());
        }
        return durations;
    }

    // Get total duration of the trajectory
    inline double getTotalDuration() const {
        double totalDuration = 0.0;
        for (int i = 0; i < getPieceNum(); i++) {
            totalDuration += pieces[i].getDuration();
        }
        return totalDuration;
    }

    inline double getTotalLength() {
        if (had_length_) {
            return length_;
        }
        length_ = 0;
        had_length_ = true;
        int num_seg = getPieceNum();
        for (size_t t = 0; t < num_seg; t++) {
            length_ += pieces[t].getLength();
        }
        return length_;

    }

    // Reload the operator[] to reach the i-th piece
    inline const Piece &operator[](int i) const {
        return pieces[i];
    }

    inline Piece &operator[](int i) {
        return pieces[i];
    }

    inline void clear(void) {
        pieces.clear();
    }

    inline Pieces::const_iterator begin() const {
        return pieces.begin();
    }

    inline Pieces::const_iterator end() const {
        return pieces.end();
    }

    // Put another piece at the tail of this trajectory
    inline void emplace_back(const Piece &piece) {
        pieces.emplace_back(piece);
        return;
    }

    // Two corresponding constructors of Piece both are supported here
    template<typename ArgTypeL, typename ArgTypeR>
    inline void emplace_back(const ArgTypeL &argL, const ArgTypeR &argR) {
        pieces.emplace_back(argL, argR);
        return;
    }

    // Append another Trajectory at the tail of this trajectory
    inline void append(const Trajectory &traj) {
        pieces.insert(pieces.end(), traj.begin(), traj.end());
        return;
    }

    // Find the piece at which the time t is located
    // The index is returned and the offset in t is removed
    inline int locatePieceIdx(double &t) const {
        int idx;
        double dur;
        for (idx = 0;
             idx < getPieceNum() &&
             t > (dur = pieces[idx].getDuration());
             idx++) {
            t -= dur;
        }
        if (idx == getPieceNum()) {
            idx--;
            t += pieces[idx].getDuration();
        }
        return idx;
    }

    // Get the position at time t of the trajectory
    inline Eigen::Vector3d getPos(double t) const {
        if (t == -1) {
            return pieces[getPieceNum() - 1].getPos(-1);
        }
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getPos(t);
    }

    // Get the velocity at time t of the trajectory
    inline Eigen::Vector3d getVel(double t) const {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getVel(t);
    }

    // Get the acceleration at time t of the trajectory
    inline Eigen::Vector3d getAcc(double t) const {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getAcc(t);
    }

    // Get the acceleration at time t of the trajectory
    inline Eigen::Vector3d getJerk(double t) const {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getJerk(t);
    }

    // Get the acceleration at time t of the trajectory
    inline Eigen::Vector3d getSnap(double t) const {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getSnap(t);
    }

    inline std::vector<Eigen::Vector3d> getUniformWaypoints(int num_per_piece = 2) {
        int num_seg = getPieceNum();
        int total_num = num_per_piece * num_seg + num_seg - 1;
        double dur = getTotalDuration();
        double seg_t = dur / total_num;
        double eval_t = seg_t;
        vector<Vec3> waypts;
        while (eval_t < dur) {
            waypts.push_back(getPos(eval_t));
            eval_t += seg_t;
            if (dur - eval_t < seg_t) {
                break;
            }
        }
        return waypts;
    }

    inline std::vector<Eigen::Vector3d> getWaypoints() {

        std::vector<Eigen::Vector3d> waypts;
        int num_seg = getPieceNum();
        for (size_t i = 1; i < num_seg; i++) {
            waypts.push_back(pieces[i].getPos(0));
        }

        /* double cur_dur, total_dur = getTotalDuration();
         for (size_t i = 0; i < num_seg; i++) {
             waypts.push_back(pieces[i].getPos(0));
             if (!pure_waypts) {
                 cur_dur = pieces[i].getDuration();
                 if (cur_dur > 0.4 * total_dur) {
                     waypts.push_back(pieces[i].getPos(cur_dur / 3));
                     waypts.push_back(pieces[i].getPos(cur_dur / 3 * 2));
                 }
             }
         }
         waypts.push_back(getPos(total_dur));*/
        return waypts;
    }


    // Get the position at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncPos(int juncIdx) const {
        if (juncIdx != getPieceNum()) {
            return pieces[juncIdx].getPos(0.0);
        } else {
            return pieces[juncIdx - 1].getPos(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the velocity at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncVel(int juncIdx) const {
        if (juncIdx != getPieceNum()) {
            return pieces[juncIdx].getVel(0.0);
        } else {
            return pieces[juncIdx - 1].getVel(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the acceleration at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncAcc(int juncIdx) const {
        if (juncIdx != getPieceNum()) {
            return pieces[juncIdx].getAcc(0.0);
        } else {
            return pieces[juncIdx - 1].getAcc(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the acceleration at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncJerk(int juncIdx) const {
        if (juncIdx != getPieceNum()) {
            return pieces[juncIdx].getJerk(0.0);
        } else {
            return pieces[juncIdx - 1].getJerk(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the max velocity rate of the trajectory
    inline double getMaxVelRate() const {
        double maxVelRate = -INFINITY;
        double tempNorm;
        for (int i = 0; i < getPieceNum(); i++) {
            tempNorm = pieces[i].getMaxVelRate();
            maxVelRate = maxVelRate < tempNorm ? tempNorm : maxVelRate;
        }
        return maxVelRate;
    }

    // Get the max velocity rate of the trajectory
    inline double getMinVelRate() const {
        double minVelRate = INFINITY;
        double tempNorm;
        for (int i = 0; i < getPieceNum(); i++) {
            tempNorm = pieces[i].getMinVelRate();
            minVelRate = minVelRate > tempNorm ? tempNorm : minVelRate;
        }
        return minVelRate;
    }

    // Get the max acceleration rate of the trajectory
    inline double getMaxAccRate() const {
        double maxAccRate = -INFINITY;
        double tempNorm;
        for (int i = 0; i < getPieceNum(); i++) {
            tempNorm = pieces[i].getMaxAccRate();
            maxAccRate = maxAccRate < tempNorm ? tempNorm : maxAccRate;
        }
        return maxAccRate;
    }

    // Check whether the velocity rate of this trajectory exceeds the threshold
    inline bool checkMaxVelRate(double maxVelRate) const {
        bool feasible = true;
        for (int i = 0; i < getPieceNum() && feasible; i++) {
            feasible = feasible && pieces[i].checkMaxVelRate(maxVelRate);
        }
        return feasible;
    }

    // Check whether the acceleration rate of this trajectory exceeds the threshold
    inline bool checkMaxAccRate(double maxAccRate) const {
        bool feasible = true;
        for (int i = 0; i < getPieceNum() && feasible; i++) {
            feasible = feasible && pieces[i].checkMaxAccRate(maxAccRate);
        }
        return feasible;
    }

    // Scale the Trajectory(t) to Trajectory(k*t)
    inline void scaleTime(double k) {
        for (int i = 0; i < getPieceNum(); i++) {
            pieces[i].scaleTime(k);
        }
    }

    inline void reserve(const int &n) {
        pieces.reserve(n);
        return;
    }

    inline void getRotation(const double &t,
                            const double &yaw,
                            const double &gAcc,
                            Eigen::Matrix3d &rotM) const {
        rotM.col(2) = getAcc(t);
        rotM(2, 2) += gAcc;
        rotM.col(2).normalize();
        rotM.col(1) = rotM.col(2).cross(Eigen::Vector3d(cos(yaw), sin(yaw), 0.0));
        rotM.col(1).normalize();
        rotM.col(0) = rotM.col(1).cross(rotM.col(2));
        return;
    }
};


// The banded system class is used for solving
// banded linear system Ax=b efficiently.
// A is an N*N band matrix with lower band width lowerBw
// and upper band width upperBw.
// Banded LU factorization has O(N) time complexity.
class BandedSystem {
public:
    // The size of A, as well as the lower/upper
    // banded width p/q are needed
    inline void create(const int &n, const int &p, const int &q) {
        // In case of re-creating before destroying
        destroy();
        N = n;
        lowerBw = p;
        upperBw = q;
        int actualSize = N * (lowerBw + upperBw + 1);
        ptrData = new double[actualSize];
        std::fill_n(ptrData, actualSize, 0.0);
        return;
    }

    inline void destroy() {
        if (ptrData != nullptr) {
            delete[] ptrData;
            ptrData = nullptr;
        }
        return;
    }

private:
    int N;
    int lowerBw;
    int upperBw;
    // Compulsory nullptr initialization here
    double *ptrData = nullptr;

public:
    // Reset the matrix to zero
    inline void reset(void) {
        std::fill_n(ptrData, N * (lowerBw + upperBw + 1), 0.0);
        return;
    }

    // The band matrix is stored as suggested in "Matrix Computation"
    inline const double &operator()(const int &i, const int &j) const {
        return ptrData[(i - j + upperBw) * N + j];
    }

    inline double &operator()(const int &i, const int &j) {
        return ptrData[(i - j + upperBw) * N + j];
    }

    // This function conducts banded LU factorization in place
    // Note that NO PIVOT is applied on the matrix "A" for efficiency!!!
    inline void factorizeLU() {
        int iM, jM;
        double cVl;
        for (int k = 0; k <= N - 2; k++) {
            iM = std::min(k + lowerBw, N - 1);
            cVl = operator()(k, k);
            for (int i = k + 1; i <= iM; i++) {
                if (operator()(i, k) != 0.0) {
                    operator()(i, k) /= cVl;
                }
            }
            jM = std::min(k + upperBw, N - 1);
            for (int j = k + 1; j <= jM; j++) {
                cVl = operator()(k, j);
                if (cVl != 0.0) {
                    for (int i = k + 1; i <= iM; i++) {
                        if (operator()(i, k) != 0.0) {
                            operator()(i, j) -= operator()(i, k) * cVl;
                        }
                    }
                }
            }
        }
        return;
    }

    // This function solves Ax=b, then stores x in b
    // The input b is required to be N*m, i.e.,
    // m vectors to be solved.
    template<typename EIGENMAT>
    inline void solve(EIGENMAT &b) const {
        int iM;
        for (int j = 0; j <= N - 1; j++) {
            iM = std::min(j + lowerBw, N - 1);
            for (int i = j + 1; i <= iM; i++) {
                if (operator()(i, j) != 0.0) {
                    b.row(i) -= operator()(i, j) * b.row(j);
                }
            }
        }
        for (int j = N - 1; j >= 0; j--) {
            b.row(j) /= operator()(j, j);
            iM = std::max(0, j - upperBw);
            for (int i = iM; i <= j - 1; i++) {
                if (operator()(i, j) != 0.0) {
                    b.row(i) -= operator()(i, j) * b.row(j);
                }
            }
        }
        return;
    }

    // This function solves ATx=b, then stores x in b
    // The input b is required to be N*m, i.e.,
    // m vectors to be solved.
    template<typename EIGENMAT>
    inline void solveAdj(EIGENMAT &b) const {
        int iM;
        for (int j = 0; j <= N - 1; j++) {
            b.row(j) /= operator()(j, j);
            iM = std::min(j + upperBw, N - 1);
            for (int i = j + 1; i <= iM; i++) {
                if (operator()(j, i) != 0.0) {
                    b.row(i) -= operator()(j, i) * b.row(j);
                }
            }
        }
        for (int j = N - 1; j >= 0; j--) {
            iM = std::max(0, j - lowerBw);
            for (int i = iM; i <= j - 1; i++) {
                if (operator()(j, i) != 0.0) {
                    b.row(i) -= operator()(j, i) * b.row(j);
                }
            }
        }
        return;
    }
};

#endif
