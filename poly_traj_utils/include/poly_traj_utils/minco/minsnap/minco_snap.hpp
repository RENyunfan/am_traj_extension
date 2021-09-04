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

#ifndef QNMOPT_HPP
#define QNMOPT_HPP

#include "root_finder.hpp"
#include "lbfgs.hpp"
#include "geoutils.hpp"

#include <Eigen/Eigen>

#include <iostream>
#include <cmath>
#include <cfloat>
#include <vector>
#include "fmt/color.h"
#include "../../poly_data_structure.hpp"




class MINCO_S4 {
public:
    MINCO_S4() = default;

    ~MINCO_S4() { A.destroy(); }

private:
    int N;
    Eigen::Matrix<double, 3, 4> headPVAJ;
    Eigen::Matrix<double, 3, 4> tailPVAJ;
    Eigen::VectorXd T1;
    BandedSystem A;
    Eigen::MatrixX3d b;

    // Temp variables
    Eigen::VectorXd T2;
    Eigen::VectorXd T3;
    Eigen::VectorXd T4;
    Eigen::VectorXd T5;
    Eigen::VectorXd T6;
    Eigen::VectorXd T7;
    Eigen::MatrixX3d gdC;
    double smoothEps;

private:
    inline void addGradJbyT(Eigen::VectorXd &gdT) const {
        for (int i = 0; i < N; i++) {
            gdT(i) += 576.0 * b.row(8 * i + 4).squaredNorm() +
                      5760.0 * b.row(8 * i + 4).dot(b.row(8 * i + 5)) * T1(i) +
                      14400.0 * b.row(8 * i + 5).squaredNorm() * T2(i) +
                      17280.0 * b.row(8 * i + 4).dot(b.row(8 * i + 6)) * T2(i) +
                      86400.0 * b.row(8 * i + 5).dot(b.row(8 * i + 6)) * T3(i) +
                      40320.0 * b.row(8 * i + 4).dot(b.row(8 * i + 7)) * T3(i) +
                      129600.0 * b.row(8 * i + 6).squaredNorm() * T4(i) +
                      201600.0 * b.row(8 * i + 5).dot(b.row(8 * i + 7)) * T4(i) +
                      604800.0 * b.row(8 * i + 6).dot(b.row(8 * i + 7)) * T5(i) +
                      705600.0 * b.row(8 * i + 7).squaredNorm() * T6(i);
        }
        return;
    }

    inline void addGradJbyC(Eigen::MatrixX3d &gdC) const {
        for (int i = 0; i < N; i++) {
            gdC.row(8 * i + 7) += 10080.0 * b.row(8 * i + 4) * T4(i) +
                                  40320.0 * b.row(8 * i + 5) * T5(i) +
                                  100800.0 * b.row(8 * i + 6) * T6(i) +
                                  201600.0 * b.row(8 * i + 7) * T7(i);

            gdC.row(8 * i + 6) += 5760.0 * b.row(8 * i + 4) * T3(i) +
                                  21600.0 * b.row(8 * i + 5) * T4(i) +
                                  51840.0 * b.row(8 * i + 6) * T5(i) +
                                  100800.0 * b.row(8 * i + 7) * T6(i);

            gdC.row(8 * i + 5) += 2880.0 * b.row(8 * i + 4) * T2(i) +
                                  9600.0 * b.row(8 * i + 5) * T3(i) +
                                  21600.0 * b.row(8 * i + 6) * T4(i) +
                                  40320.0 * b.row(8 * i + 7) * T5(i);

            gdC.row(8 * i + 4) += 1152.0 * b.row(8 * i + 4) * T1(i) +
                                  2880.0 * b.row(8 * i + 5) * T2(i) +
                                  5760.0 * b.row(8 * i + 6) * T3(i) +
                                  10080.0 * b.row(8 * i + 7) * T4(i);
        }
        return;
    }

    inline void solveAdjGradC(Eigen::MatrixX3d &gdC) const {
        A.solveAdj(gdC);
        return;
    }

    inline void addPropCtoT(const Eigen::MatrixX3d &adjGdC, Eigen::VectorXd &gdT) const {

        Eigen::Matrix<double, 8, 3> B1;
        Eigen::Matrix<double, 4, 3> B2;

        for (int i = 0; i < N - 1; i++) {
            // negative velocity
            B1.row(3) = -(b.row(i * 8 + 1) +
                          2.0 * T1(i) * b.row(i * 8 + 2) +
                          3.0 * T2(i) * b.row(i * 8 + 3) +
                          4.0 * T3(i) * b.row(i * 8 + 4) +
                          5.0 * T4(i) * b.row(i * 8 + 5) +
                          6.0 * T5(i) * b.row(i * 8 + 6) +
                          7.0 * T6(i) * b.row(i * 8 + 7));
            B1.row(4) = B1.row(3);

            // negative acceleration
            B1.row(5) = -(2.0 * b.row(i * 8 + 2) +
                          6.0 * T1(i) * b.row(i * 8 + 3) +
                          12.0 * T2(i) * b.row(i * 8 + 4) +
                          20.0 * T3(i) * b.row(i * 8 + 5) +
                          30.0 * T4(i) * b.row(i * 8 + 6) +
                          42.0 * T5(i) * b.row(i * 8 + 7));

            // negative jerk
            B1.row(6) = -(6.0 * b.row(i * 8 + 3) +
                          24.0 * T1(i) * b.row(i * 8 + 4) +
                          60.0 * T2(i) * b.row(i * 8 + 5) +
                          120.0 * T3(i) * b.row(i * 8 + 6) +
                          210.0 * T4(i) * b.row(i * 8 + 7));

            // negative snap
            B1.row(7) = -(24.0 * b.row(i * 8 + 4) +
                          120.0 * T1(i) * b.row(i * 8 + 5) +
                          360.0 * T2(i) * b.row(i * 8 + 6) +
                          840.0 * T3(i) * b.row(i * 8 + 7));

            // negative crackle
            B1.row(0) = -(120.0 * b.row(i * 8 + 5) +
                          720.0 * T1(i) * b.row(i * 8 + 6) +
                          2520.0 * T2(i) * b.row(i * 8 + 7));

            // negative d_crackle
            B1.row(1) = -(720.0 * b.row(i * 8 + 6) +
                          5040.0 * T1(i) * b.row(i * 8 + 7));

            // negative dd_crackle
            B1.row(2) = -5040.0 * b.row(i * 8 + 7);

            gdT(i) += B1.cwiseProduct(adjGdC.block<8, 3>(8 * i + 4, 0)).sum();
        }

        // negative velocity
        B2.row(0) = -(b.row(8 * N - 7) +
                      2.0 * T1(N - 1) * b.row(8 * N - 6) +
                      3.0 * T2(N - 1) * b.row(8 * N - 5) +
                      4.0 * T3(N - 1) * b.row(8 * N - 4) +
                      5.0 * T4(N - 1) * b.row(8 * N - 3) +
                      6.0 * T5(N - 1) * b.row(8 * N - 2) +
                      7.0 * T6(N - 1) * b.row(8 * N - 1));

        // negative acceleration
        B2.row(1) = -(2.0 * b.row(8 * N - 6) +
                      6.0 * T1(N - 1) * b.row(8 * N - 5) +
                      12.0 * T2(N - 1) * b.row(8 * N - 4) +
                      20.0 * T3(N - 1) * b.row(8 * N - 3) +
                      30.0 * T4(N - 1) * b.row(8 * N - 2) +
                      42.0 * T5(N - 1) * b.row(8 * N - 1));

        // negative jerk
        B2.row(2) = -(6.0 * b.row(8 * N - 5) +
                      24.0 * T1(N - 1) * b.row(8 * N - 4) +
                      60.0 * T2(N - 1) * b.row(8 * N - 3) +
                      120.0 * T3(N - 1) * b.row(8 * N - 2) +
                      210.0 * T4(N - 1) * b.row(8 * N - 1));

        // negative snap
        B2.row(3) = -(24.0 * b.row(8 * N - 4) +
                      120.0 * T1(N - 1) * b.row(8 * N - 3) +
                      360.0 * T2(N - 1) * b.row(8 * N - 2) +
                      840.0 * T3(N - 1) * b.row(8 * N - 1));

        gdT(N - 1) += B2.cwiseProduct(adjGdC.block<4, 3>(8 * N - 4, 0)).sum();

        return;
    }

    inline void addPropCtoP(const Eigen::MatrixX3d &adjGdC, Eigen::Matrix3Xd &gdInP) const {
        for (int i = 0; i < N - 1; i++) {
            gdInP.col(i) += adjGdC.row(8 * i + 7).transpose();
        }
        return;
    }

    inline void positiveSmoothedL1(const double &x, double &f, double &df) const {
        const double pe = smoothEps;
        const double half = 0.5 * pe;
        const double f3c = 1.0 / (pe * pe);
        const double f4c = -0.5 * f3c / pe;
        const double d2c = 3.0 * f3c;
        const double d3c = 4.0 * f4c;

        if (x < pe) {
            f = (f4c * x + f3c) * x * x * x;
            df = (d3c * x + d2c) * x * x;
        } else {
            f = x - half;
            df = 1.0;
        }

        return;
    }

    inline void addTimeIntPenalty(const Eigen::VectorXi cons,       /// cons = 8 是用来拆分轨迹作离散化的
                                  const double vmax,    /// 最大速度限制
                                  const double amax,    /// 最大加速度读限制
                                  const Eigen::Vector3d ci, /// 权重系数
                                  double &cost,         /// 代价函数
                                  Eigen::VectorXd &gdT, /// t的梯度向量
                                  Eigen::MatrixX3d &gdC) const {


        double pena = 0.0;
        const double vmaxSqr = vmax * vmax;
        const double amaxSqr = amax * amax;

        Eigen::Vector3d pos, vel, acc, jer;
        double step, alpha;
        double s1, s2, s3, s4, s5, s6, s7;
        Eigen::Matrix<double, 8, 1> beta0, beta1, beta2, beta3;
        Eigen::Vector3d outerNormal;
        int K;
        double violaPos, violaVel, violaAcc;
        double violaPosPenaD, violaVelPenaD, violaAccPenaD;
        double violaPosPena, violaVelPena, violaAccPena;
        Eigen::Matrix<double, 8, 3> gradViolaVc, gradViolaAc;
        double gradViolaVt, gradViolaAt;
        double omg;

        int innerLoop, idx;
        for (int i = 0; i < N; i++) {
            // c是多项式系数
            const auto &c = b.block<8, 3>(i * 8, 0);
            /* cons
             *  大小为轨迹段数
             *  内容为一组与优化参数
             *  代码中设值为8
             *  */
            step = T1(i) / cons(i);
            s1 = 0.0;
            innerLoop = cons(i) + 1;

            // innerLoop的次数是固定的8+1=9 应该是用于将每段轨迹采样9次
            for (int j = 0; j < innerLoop; j++) {
                // s1是拆分八段后的轨迹，也就是单段轨迹按八分之sample来计算cost。
                s2 = s1 * s1;
                s3 = s2 * s1;
                s4 = s2 * s2;
                s5 = s4 * s1;
                s6 = s4 * s2;
                s7 = s4 * s3;
                beta0(0) = 1.0, beta0(1) = s1, beta0(2) = s2, beta0(3) = s3, beta0(4) = s4, beta0(5) = s5, beta0(
                        6) = s6, beta0(7) = s7;
                beta1(0) = 0.0, beta1(1) = 1.0, beta1(2) = 2.0 * s1, beta1(3) = 3.0 * s2, beta1(4) =
                        4.0 * s3, beta1(
                        5) = 5.0 * s4, beta1(6) = 6.0 * s5, beta1(7) = 7.0 * s6;
                beta2(0) = 0.0, beta2(1) = 0.0, beta2(2) = 2.0, beta2(3) = 6.0 * s1, beta2(4) = 12.0 * s2, beta2(
                        5) =
                        20.0 * s3, beta2(6) = 30.0 * s4, beta2(7) = 42.0 * s5;
                beta3(0) = 0.0, beta3(1) = 0.0, beta3(2) = 0.0, beta3(3) = 6.0, beta3(4) = 24.0 * s1, beta3(5) =
                        60.0 * s2, beta3(6) = 120.0 * s3, beta3(7) = 210.0 * s4;
                alpha = 1.0 / cons(i) * j;
                pos = c.transpose() * beta0;
                vel = c.transpose() * beta1;
                acc = c.transpose() * beta2;
                jer = c.transpose() * beta3;
                violaVel = vel.squaredNorm() - vmaxSqr;
                violaAcc = acc.squaredNorm() - amaxSqr;

                omg = (j == 0 || j == innerLoop - 1) ? 0.5 : 1.0;


                if (violaVel > 0.0) {
                    positiveSmoothedL1(violaVel, violaVelPena, violaVelPenaD);
                    gradViolaVc = 2.0 * beta1 * vel.transpose();
                    gradViolaVt = 2.0 * alpha * vel.transpose() * acc;
                    gdC.block<8, 3>(i * 8, 0) += omg * step * ci(1) * violaVelPenaD * gradViolaVc;
                    gdT(i) += omg * (ci(1) * violaVelPenaD * gradViolaVt * step +
                                     ci(1) * violaVelPena / cons(i));
                    pena += omg * step * ci(1) * violaVelPena;
                }

                if (violaAcc > 0.0) {
                    positiveSmoothedL1(violaAcc, violaAccPena, violaAccPenaD);
                    gradViolaAc = 2.0 * beta2 * acc.transpose();
                    gradViolaAt = 2.0 * alpha * acc.transpose() * jer;
                    gdC.block<8, 3>(i * 8, 0) += omg * step * ci(2) * violaAccPenaD * gradViolaAc;
                    gdT(i) += omg * (ci(2) * violaAccPenaD * gradViolaAt * step +
                                     ci(2) * violaAccPena / cons(i));
                    pena += omg * step * ci(2) * violaAccPena;
                }

                s1 += step;
            }
        }

        cost += pena;
        return;
    }

public:
    inline void reset(const StatePVAJ &headState,
                      const StatePVAJ &tailState,
                      const double &smoEps,
                      const int &pieceNum) {
        // 初始话轨迹段数
        N = pieceNum;
        // 初始化边界条件
        headPVAJ = headState;
        tailPVAJ = tailState;
        // 初始化时间分配向量
        T1.resize(N);
        // 初始化轨迹生成LU分解矩阵
        A.create(8 * N, 8, 8);
        b.resize(8 * N, 3);
        // 初始化多项式系数梯度
        gdC.resize(8 * N, 3);
        // 初始化L1 smooth的参数
        smoothEps = smoEps;
        return;
    }

    inline void generate(const Eigen::MatrixXd &inPs,
                         const Eigen::VectorXd &ts) {
        T1 = ts;
        T2 = T1.cwiseProduct(T1);
        T3 = T2.cwiseProduct(T1);
        T4 = T2.cwiseProduct(T2);
        T5 = T4.cwiseProduct(T1);
        T6 = T4.cwiseProduct(T2);
        T7 = T4.cwiseProduct(T3);

        A.reset();
        b.setZero();

        A(0, 0) = 1.0;
        A(1, 1) = 1.0;
        A(2, 2) = 2.0;
        A(3, 3) = 6.0;
        b.row(0) = headPVAJ.col(0).transpose();
        b.row(1) = headPVAJ.col(1).transpose();
        b.row(2) = headPVAJ.col(2).transpose();
        b.row(3) = headPVAJ.col(3).transpose();

        for (int i = 0; i < N - 1; i++) {
            A(8 * i + 4, 8 * i + 4) = 24.0;
            A(8 * i + 4, 8 * i + 5) = 120.0 * T1(i);
            A(8 * i + 4, 8 * i + 6) = 360.0 * T2(i);
            A(8 * i + 4, 8 * i + 7) = 840.0 * T3(i);
            A(8 * i + 4, 8 * i + 12) = -24.0;

            A(8 * i + 5, 8 * i + 5) = 120.0;
            A(8 * i + 5, 8 * i + 6) = 720.0 * T1(i);
            A(8 * i + 5, 8 * i + 7) = 2520.0 * T2(i);
            A(8 * i + 5, 8 * i + 13) = -120.0;

            A(8 * i + 6, 8 * i + 6) = 720.0;
            A(8 * i + 6, 8 * i + 7) = 5040.0 * T1(i);
            A(8 * i + 6, 8 * i + 14) = -720.0;

            A(8 * i + 7, 8 * i) = 1.0;
            A(8 * i + 7, 8 * i + 1) = T1(i);
            A(8 * i + 7, 8 * i + 2) = T2(i);
            A(8 * i + 7, 8 * i + 3) = T3(i);
            A(8 * i + 7, 8 * i + 4) = T4(i);
            A(8 * i + 7, 8 * i + 5) = T5(i);
            A(8 * i + 7, 8 * i + 6) = T6(i);
            A(8 * i + 7, 8 * i + 7) = T7(i);

            A(8 * i + 8, 8 * i) = 1.0;
            A(8 * i + 8, 8 * i + 1) = T1(i);
            A(8 * i + 8, 8 * i + 2) = T2(i);
            A(8 * i + 8, 8 * i + 3) = T3(i);
            A(8 * i + 8, 8 * i + 4) = T4(i);
            A(8 * i + 8, 8 * i + 5) = T5(i);
            A(8 * i + 8, 8 * i + 6) = T6(i);
            A(8 * i + 8, 8 * i + 7) = T7(i);
            A(8 * i + 8, 8 * i + 8) = -1.0;

            A(8 * i + 9, 8 * i + 1) = 1.0;
            A(8 * i + 9, 8 * i + 2) = 2.0 * T1(i);
            A(8 * i + 9, 8 * i + 3) = 3.0 * T2(i);
            A(8 * i + 9, 8 * i + 4) = 4.0 * T3(i);
            A(8 * i + 9, 8 * i + 5) = 5.0 * T4(i);
            A(8 * i + 9, 8 * i + 6) = 6.0 * T5(i);
            A(8 * i + 9, 8 * i + 7) = 7.0 * T6(i);
            A(8 * i + 9, 8 * i + 9) = -1.0;

            A(8 * i + 10, 8 * i + 2) = 2.0;
            A(8 * i + 10, 8 * i + 3) = 6.0 * T1(i);
            A(8 * i + 10, 8 * i + 4) = 12.0 * T2(i);
            A(8 * i + 10, 8 * i + 5) = 20.0 * T3(i);
            A(8 * i + 10, 8 * i + 6) = 30.0 * T4(i);
            A(8 * i + 10, 8 * i + 7) = 42.0 * T5(i);
            A(8 * i + 10, 8 * i + 10) = -2.0;

            A(8 * i + 11, 8 * i + 3) = 6.0;
            A(8 * i + 11, 8 * i + 4) = 24.0 * T1(i);
            A(8 * i + 11, 8 * i + 5) = 60.0 * T2(i);
            A(8 * i + 11, 8 * i + 6) = 120.0 * T3(i);
            A(8 * i + 11, 8 * i + 7) = 210.0 * T4(i);
            A(8 * i + 11, 8 * i + 11) = -6.0;

            b.row(8 * i + 7) = inPs.col(i).transpose();
        }

        A(8 * N - 4, 8 * N - 8) = 1.0;
        A(8 * N - 4, 8 * N - 7) = T1(N - 1);
        A(8 * N - 4, 8 * N - 6) = T2(N - 1);
        A(8 * N - 4, 8 * N - 5) = T3(N - 1);
        A(8 * N - 4, 8 * N - 4) = T4(N - 1);
        A(8 * N - 4, 8 * N - 3) = T5(N - 1);
        A(8 * N - 4, 8 * N - 2) = T6(N - 1);
        A(8 * N - 4, 8 * N - 1) = T7(N - 1);

        A(8 * N - 3, 8 * N - 7) = 1.0;
        A(8 * N - 3, 8 * N - 6) = 2.0 * T1(N - 1);
        A(8 * N - 3, 8 * N - 5) = 3.0 * T2(N - 1);
        A(8 * N - 3, 8 * N - 4) = 4.0 * T3(N - 1);
        A(8 * N - 3, 8 * N - 3) = 5.0 * T4(N - 1);
        A(8 * N - 3, 8 * N - 2) = 6.0 * T5(N - 1);
        A(8 * N - 3, 8 * N - 1) = 7.0 * T6(N - 1);

        A(8 * N - 2, 8 * N - 6) = 2.0;
        A(8 * N - 2, 8 * N - 5) = 6.0 * T1(N - 1);
        A(8 * N - 2, 8 * N - 4) = 12.0 * T2(N - 1);
        A(8 * N - 2, 8 * N - 3) = 20.0 * T3(N - 1);
        A(8 * N - 2, 8 * N - 2) = 30.0 * T4(N - 1);
        A(8 * N - 2, 8 * N - 1) = 42.0 * T5(N - 1);

        A(8 * N - 1, 8 * N - 5) = 6.0;
        A(8 * N - 1, 8 * N - 4) = 24.0 * T1(N - 1);
        A(8 * N - 1, 8 * N - 3) = 60.0 * T2(N - 1);
        A(8 * N - 1, 8 * N - 2) = 120.0 * T3(N - 1);
        A(8 * N - 1, 8 * N - 1) = 210.0 * T4(N - 1);

        b.row(8 * N - 4) = tailPVAJ.col(0).transpose();
        b.row(8 * N - 3) = tailPVAJ.col(1).transpose();
        b.row(8 * N - 2) = tailPVAJ.col(2).transpose();
        b.row(8 * N - 1) = tailPVAJ.col(3).transpose();

        A.factorizeLU();
        A.solve(b);

        return;
    }

    inline double getTrajSnapCost() const {
        double objective = 0.0;
        for (int i = 0; i < N; i++) {
            objective += 576.0 * b.row(8 * i + 4).squaredNorm() * T1(i) +
                         2880.0 * b.row(8 * i + 4).dot(b.row(8 * i + 5)) * T2(i) +
                         4800.0 * b.row(8 * i + 5).squaredNorm() * T3(i) +
                         5760.0 * b.row(8 * i + 4).dot(b.row(8 * i + 6)) * T3(i) +
                         21600.0 * b.row(8 * i + 5).dot(b.row(8 * i + 6)) * T4(i) +
                         10080.0 * b.row(8 * i + 4).dot(b.row(8 * i + 7)) * T4(i) +
                         25920.0 * b.row(8 * i + 6).squaredNorm() * T5(i) +
                         40320.0 * b.row(8 * i + 5).dot(b.row(8 * i + 7)) * T5(i) +
                         100800.0 * b.row(8 * i + 6).dot(b.row(8 * i + 7)) * T6(i) +
                         100800.0 * b.row(8 * i + 7).squaredNorm() * T7(i);
        }
        return objective;
    }

    inline void evalTrajCostGrad(const Eigen::VectorXi &cons,
                                 const double &vmax,
                                 const double &amax,
                                 const Eigen::Vector3d &ci,
                                 double &cost,
                                 Eigen::VectorXd &gdT,
                                 Eigen::Matrix3Xd &gdInPs) {
        gdT.setZero();
        gdC.setZero();
        gdInPs.setZero();
        cost = 0;
        // 首先计算轨迹的snapcost
        cost = getTrajSnapCost();
        // 计算snap代价对于时间梯度
        addGradJbyT(gdT);
        // 计算snap代价对于多项式系数梯度
        addGradJbyC(gdC);
        // 计算积分代价和积分梯度                              输入gdT和gdC
        addTimeIntPenalty(cons, vmax, amax, ci, cost, gdT, gdC);
        solveAdjGradC(gdC);
        addPropCtoT(gdC, gdT);
        //        Eigen::Matrix3Xd gdInPsfull;
        //        int id = 0;
        //        for(int i = 0 ; i < gdInPsfull.cols()){
        //            if(i%2){
        //                gdInPs(id) = gdInPsfull(i);
        //                id++;
        //            }
        //        }
        addPropCtoP(gdC, gdInPs);

    }

    inline Trajectory getTraj(void) const {
        Trajectory traj;
        traj.reserve(N);
        for (int i = 0; i < N; i++) {
            traj.emplace_back(T1(i), b.block<8, 3>(8 * i, 0).transpose().rowwise().reverse());
        }
        return traj;
    }
};

class QNMOPT {
private:
    // Use C2 or Cinf diffeo
    bool c2dfm;

    // Use retraction for hyperball
    bool retraction;

    // Use soft time or not
    bool softT;

    // Weight for time regularization term
    double rho;

    // Fixed total time
    double sumT;

    //Minimum Jerk Optimizer
    MINCO_S4 minco;

    // Temp variables for problem solving
    StatePVAJ iState;
    StatePVAJ fState;
    vector<Vec3> wayPts;

    // Each col of cfgHs denotes a facet (outter_normal^T,point^T)^T
    std::vector<Eigen::Matrix3Xd> cfgVs;
    std::vector<Eigen::Matrix<double, 6, -1>> cfgHs;
    Eigen::Matrix3Xd gdInPs;

    // Piece num for each polytope
    Eigen::VectorXi intervals;
    // Assignment vector for point in V-polytope
    Eigen::VectorXi idxVs;
    // Assignment vector for piece in H-polytope
    Eigen::VectorXi idxHs;

    int coarseN;
    int fineN;
    int dimFreeT;
    int dimFreeP;
    Eigen::VectorXd coarseT;
    Eigen::VectorXd fineT;
    Eigen::Matrix3Xd innerP;
    Eigen::Matrix3Xd freePts;
    // Params for constraints
    Eigen::VectorXi cons;
    double smoothEps;
    Eigen::Vector3d chi;
    double vmax;
    double amax;

    // L-BFGS Solver Parameters
    lbfgs::lbfgs_parameter_t lbfgs_params;

private:
    template<typename EIGENVEC>
    static inline void forwardT(const EIGENVEC &t,
                                Eigen::VectorXd &vecT,
                                bool soft,
                                const double &sT,
                                bool c2) {
        if (soft) {
            if (c2) {
                int M = vecT.size();
                for (int i = 0; i < M; i++) {
                    vecT(i) = t(i) > 0.0
                              ? ((0.5 * t(i) + 1.0) * t(i) + 1.0)
                              : 1.0 / ((0.5 * t(i) - 1.0) * t(i) + 1.0);
                }
            } else {
                vecT = t.array().exp();
            }
        } else {
            if (c2) {
                int Ms1 = t.size();
                for (int i = 0; i < Ms1; i++) {
                    vecT(i) = t(i) > 0.0
                              ? ((0.5 * t(i) + 1.0) * t(i) + 1.0)
                              : 1.0 / ((0.5 * t(i) - 1.0) * t(i) + 1.0);
                }
                vecT(Ms1) = 0.0;
                vecT /= 1.0 + vecT.sum();
                vecT(Ms1) = 1.0 - vecT.sum();
                vecT *= sT;
            } else {
                int Ms1 = t.size();
                vecT.head(Ms1) = t.array().exp();
                vecT(Ms1) = 0.0;
                vecT /= 1.0 + vecT.sum();
                vecT(Ms1) = 1.0 - vecT.sum();
                vecT *= sT;
            }
        }
        return;
    }

    template<typename EIGENVEC>
    static inline void backwardT(const Eigen::VectorXd &vecT, ///<[in,out] 初始时间分配
                                 EIGENVEC &t,   ///<[in,out] 时间状态向量
                                 bool soft,     ///<[in] 是否使用SOFT时间
                                 bool c2      ///<[in] 是否使用c2连续
    ) {

        if (soft) {
            if (c2) {
                int M = vecT.size();
                for (int i = 0; i < M; i++) {
                    // 更新优化器的初始状态时间为 一种归一化的形式
                    t(i) = vecT(i) > 1.0
                           ? (sqrt(2.0 * vecT(i) - 1.0) - 1.0)
                           : (1.0 - sqrt(2.0 / vecT(i) - 1.0));
                }
            } else {
                t = vecT.array().log();
            }
        } else {
            if (c2) {
                int Ms1 = t.size();
                t = vecT.head(Ms1) / vecT(Ms1);
                for (int i = 0; i < Ms1; i++) {
                    t(i) = t(i) > 1.0
                           ? (sqrt(2.0 * t(i) - 1.0) - 1.0)
                           : (1.0 - sqrt(2.0 / t(i) - 1.0));
                }
            } else {
                int Ms1 = t.size();
                t = (vecT.head(Ms1) / vecT(Ms1)).array().log();
            }
        }
        return;
    }

    /// forwardP
    template<typename EIGENVEC>
    static inline void forwardP(const EIGENVEC &p,      ///<[in] 输入状态向量
                                const Eigen::VectorXi &idVs,    ///<[in] 输入顶点状态
                                const std::vector<Eigen::Matrix3Xd> &cfgPolyVs,///<[in] 输入多项式
                                Eigen::Matrix3Xd &inP,///<[in]    输出位置
                                bool retract        ///<[in]    输入是否使用retract
    ) {
        if (!retract) {
            int M = inP.cols();
            Eigen::VectorXd q;
            int j = 0, k, idx;
            for (int i = 0; i < M; i++) {
                idx = idVs(i);
                k = cfgPolyVs[idx].cols() - 1;
                q = 2.0 / (1.0 + p.segment(j, k).squaredNorm()) * p.segment(j, k);
                inP.col(i) = cfgPolyVs[idx].rightCols(k) * q.cwiseProduct(q) +
                             cfgPolyVs[idx].col(0);
                j += k;
            }
        } else {
            int M = inP.cols();
            Eigen::VectorXd q;
            double normInv;
            int j = 0, k, idx;
            for (int i = 0; i < M; i++) {
                idx = idVs(i);
                k = cfgPolyVs[idx].cols() - 1;
                normInv = 1.0 / p.segment(j, k + 1).norm();
                q = p.segment(j, k + 1).head(k) * normInv;
                inP.col(i) = cfgPolyVs[idx].rightCols(k) * q.cwiseProduct(q) +
                             cfgPolyVs[idx].col(0);
                j += k + 1;
            }
        }
        return;
    }

    static inline double objectiveNLS(void *ptrPOBs,
                                      const double *x,
                                      double *grad,
                                      const int n) {
        // 这里是backwardP中的非线性最小二乘法维度目标函数，pob是[waypoint, k个顶点]
        const Eigen::Matrix3Xd &pobs = *(Eigen::Matrix3Xd *) ptrPOBs;
        // 将n个优化变量映射到p
        Eigen::Map<const Eigen::VectorXd> p(x, n);
        // 将梯度映射到gradp
        Eigen::Map<Eigen::VectorXd> gradp(grad, n);
        //  p一开始还是空的
        double qnsqr = p.squaredNorm();
        double qnsqrp1 = qnsqr + 1.0;
        double qnsqrp1sqr = qnsqrp1 * qnsqrp1;
        Eigen::VectorXd r = 2.0 / qnsqrp1 * p;
        //delta等于r和每一顶点的内
        Eigen::Vector3d delta = pobs.rightCols(n) * r.cwiseProduct(r) +
                                pobs.col(1) - pobs.col(0);
        double cost = delta.squaredNorm();
        Eigen::Vector3d gradR3 = 2 * delta;
        // 根据公式82计算梯度
        Eigen::VectorXd gdr = pobs.rightCols(n).transpose() * gradR3;
        gdr = gdr.array() * r.array() * 2.0;
        gradp = gdr * 2.0 / qnsqrp1 -
                p * 4.0 * gdr.dot(p) / qnsqrp1sqr;

        return cost;
    }

    static inline double objectiveRetractNLS(void *ptrPOBs,
                                             const double *x,
                                             double *grad,
                                             const int n) {
        const Eigen::Matrix3Xd &pobs = *(Eigen::Matrix3Xd *) ptrPOBs;
        Eigen::Map<const Eigen::VectorXd> p(x, n);
        Eigen::Map<Eigen::VectorXd> gradp(grad, n);

        double sqrnormp = p.squaredNorm();
        double normInv = 1.0 / sqrt(sqrnormp);
        Eigen::VectorXd unitp = p * normInv;
        Eigen::VectorXd r = unitp.head(n - 1);

        Eigen::Vector3d delta = pobs.rightCols(n - 1) * r.cwiseProduct(r) +
                                pobs.col(1) - pobs.col(0);
        double cost = delta.squaredNorm();
        Eigen::Vector3d gradR3 = 2 * delta;

        Eigen::VectorXd gdr = pobs.rightCols(n - 1).transpose() * gradR3;
        gdr = gdr.array() * r.array() * 2.0;

        gradp.head(n - 1) = gdr;
        gradp(n - 1) = 0.0;
        gradp = (gradp - unitp * unitp.dot(gradp)).eval() * normInv;

        const double normViola = sqrnormp - 1.0;
        if (normViola > 0.0) {
            double c = normViola * normViola;
            const double dc = 3.0 * c;
            c *= normViola;
            cost += c;
            gradp += dc * 2.0 * p;
        }

        return cost;
    }

    template<typename EIGENVEC>
    static inline void backwardP(const Eigen::Matrix3Xd &inP, ///<[in] 初始化的waypoint
                                 const Eigen::VectorXi &idVs, ///<[in] 飞行走廊的顶点坐标
                                 const std::vector<Eigen::Matrix3Xd> &cfgPolyVs, ///<[in] 飞行走廊的顶点们
                                 EIGENVEC &p,  ///<[in] 优化的状态{目前就是顶点们}
                                 bool retract   ///<[in] 是否需要retrack
    ) {
        int M = inP.cols();
        int j = 0, k, idx;

        // Parameters for tiny nonlinear least squares
        // 进行一个小型的非线性最小二乘
        double minSqrD;
        lbfgs::lbfgs_parameter_t nls_params;
        lbfgs::lbfgs_load_default_parameters(&nls_params);
        // FLT_EPSILON 是float精度下的最小值
        nls_params.g_epsilon = FLT_EPSILON;
        nls_params.max_iterations = 128;

        // pobs目前还不知道是干啥的，也没分配内存
        Eigen::Matrix3Xd pobs;
        if (!retract) {
            for (int i = 0; i < M; i++) {
                // 获取交界处顶点id
                idx = idVs(i);
                // k是交界处的顶点个数-1
                k = cfgPolyVs[idx].cols() - 1;
                // p的 j-k是第i个走廊的优化变量
                p.segment(j, k).setConstant(1.0 / (sqrt(k + 1.0) + 1.0));
                // 初始化为3*顶点个数+1个）
                pobs.resize(3, k + 2);
                // 第一个点是需要优化的 waypoint
                pobs.col(0) = inP.col(i);
                // 剩下的是k个顶点
                pobs.rightCols(k + 1) = cfgPolyVs[idx];
                lbfgs::lbfgs_optimize(k,                // 优化变量有k个，等于顶点个数
                                      p.data() + j,     // 初始话的优化变量就是已经设置好的
                                      &minSqrD,
                                      &QNMOPT::objectiveNLS,
                                      nullptr,
                                      nullptr,
                                      &pobs,
                                      &nls_params);

                j += k;
            }
        } else {
            // pass
            for (int i = 0; i < M; i++) {
                idx = idVs(i);
                k = cfgPolyVs[idx].cols() - 1;
                p.segment(j, k + 1).setConstant(1.0 / sqrt(k + 1.0));
                pobs.resize(3, k + 2);
                pobs.col(0) = inP.col(i);
                pobs.rightCols(k + 1) = cfgPolyVs[idx];
                lbfgs::lbfgs_optimize(k + 1,
                                      p.data() + j,
                                      &minSqrD,
                                      &QNMOPT::objectiveRetractNLS,
                                      nullptr,
                                      nullptr,
                                      &pobs,
                                      &nls_params);

                j += k + 1;
            }
        }

        return;
    }

    static inline void addLayerTGrad(const Eigen::VectorXd &t,
                                     Eigen::VectorXd &gradT,
                                     bool soft,
                                     const double &sT,
                                     bool c2) {
        if (soft) {
            if (c2) {
                int M = t.size();
                double denSqrt;
                for (int i = 0; i < M; i++) {
                    if (t(i) > 0) {
                        gradT(i) *= t(i) + 1.0;
                    } else {
                        denSqrt = (0.5 * t(i) - 1.0) * t(i) + 1.0;
                        gradT(i) *= (1.0 - t(i)) / (denSqrt * denSqrt);
                    }
                }
            } else {
                int M = t.size();
                gradT.head(M).array() *= t.array().exp();
            }
        } else {
            if (c2) {
                int Ms1 = t.size();
                Eigen::VectorXd gFree = sT * gradT.head(Ms1);
                double gTail = sT * gradT(Ms1);
                Eigen::VectorXd dExpTau(Ms1);
                double expTauSum = 0.0, gFreeDotExpTau = 0.0;
                double denSqrt, expTau;
                for (int i = 0; i < Ms1; i++) {
                    if (t(i) > 0) {
                        expTau = (0.5 * t(i) + 1.0) * t(i) + 1.0;
                        dExpTau(i) = t(i) + 1.0;
                        expTauSum += expTau;
                        gFreeDotExpTau += expTau * gFree(i);
                    } else {
                        denSqrt = (0.5 * t(i) - 1.0) * t(i) + 1.0;
                        expTau = 1.0 / denSqrt;
                        dExpTau(i) = (1.0 - t(i)) / (denSqrt * denSqrt);
                        expTauSum += expTau;
                        gFreeDotExpTau += expTau * gFree(i);
                    }
                }
                denSqrt = expTauSum + 1.0;
                gradT.head(Ms1) = (gFree.array() - gTail) * dExpTau.array() / denSqrt -
                                  (gFreeDotExpTau - gTail * expTauSum) * dExpTau.array() / (denSqrt * denSqrt);
                gradT(Ms1) = 0.0;
            } else {
                int Ms1 = t.size();
                Eigen::VectorXd gFree = sT * gradT.head(Ms1);
                double gTail = sT * gradT(Ms1);
                Eigen::VectorXd expTau = t.array().exp();
                double expTauSum = expTau.sum();
                double denom = expTauSum + 1.0;
                gradT.head(Ms1) = (gFree.array() - gTail) * expTau.array() / denom -
                                  (gFree.dot(expTau) - gTail * expTauSum) * expTau.array() / (denom * denom);
                gradT(Ms1) = 0.0;
            }
        }
        return;
    }

    template<typename EIGENVEC>
    static inline void addLayerPGrad(const Eigen::VectorXd &p,
                                     const Eigen::VectorXi &idVs,
                                     const std::vector<Eigen::Matrix3Xd> &cfgPolyVs,
                                     const Eigen::Matrix3Xd &gradInPs,
                                     EIGENVEC &grad,
                                     bool retract) {
        int M = gradInPs.cols();

        if (!retract) {
            int j = 0, k, idx;
            double qnsqr, qnsqrp1, qnsqrp1sqr;
            Eigen::VectorXd q, r, gdr;
            for (int i = 0; i < M; i++) {
                idx = idVs(i);
                k = cfgPolyVs[idx].cols() - 1;

                q = p.segment(j, k);
                qnsqr = q.squaredNorm();
                qnsqrp1 = qnsqr + 1.0;
                qnsqrp1sqr = qnsqrp1 * qnsqrp1;
                r = 2.0 / qnsqrp1 * q;
                gdr = cfgPolyVs[idx].rightCols(k).transpose() * gradInPs.col(i);
                gdr = gdr.array() * r.array() * 2.0;

                grad.segment(j, k) = gdr * 2.0 / qnsqrp1 -
                                     q * 4.0 * gdr.dot(q) / qnsqrp1sqr;

                j += k;
            }
        } else {
            int j = 0, k, idx;
            double normInv;
            Eigen::VectorXd q, r, gdr, gradq, unitq;
            for (int i = 0; i < M; i++) {
                idx = idVs(i);
                k = cfgPolyVs[idx].cols() - 1;

                q = p.segment(j, k + 1);
                normInv = 1.0 / q.norm();
                unitq = q * normInv;
                r = unitq.head(k);
                gdr = cfgPolyVs[idx].rightCols(k).transpose() * gradInPs.col(i);
                gdr = gdr.array() * r.array() * 2.0;

                gradq.resize(k + 1);
                gradq.head(k) = gdr;
                gradq(k) = 0.0;
                grad.segment(j, k + 1) = (gradq - unitq * unitq.dot(gradq)) * normInv;

                j += k + 1;
            }
        }

        return;
    }

    template<typename EIGENVEC>
    static inline void restrictNorm(const Eigen::VectorXd &p,
                                    const Eigen::VectorXi &idVs,
                                    const std::vector<Eigen::Matrix3Xd> &cfgPolyVs,
                                    double &cost,
                                    EIGENVEC &grad) {
        int M = idVs.size();

        int j = 0, k, idx;
        double sqrnormq, normViola, c, dc;
        Eigen::VectorXd q;
        for (int i = 0; i < M; i++) {
            idx = idVs(i);
            k = cfgPolyVs[idx].cols() - 1;

            q = p.segment(j, k + 1);
            sqrnormq = q.squaredNorm();
            normViola = sqrnormq - 1.0;
            if (normViola > 0.0) {
                c = normViola * normViola;
                dc = 3.0 * c;
                c *= normViola;
                cost += c;
                grad.segment(j, k + 1) += dc * 2.0 * q;
            }

            j += k + 1;
        }

        return;
    }

    static inline void splitToFineT(const Eigen::VectorXd &cT,
                                    const Eigen::VectorXi &intervs,
                                    Eigen::VectorXd &fT) {
        int M = intervs.size();
        int offset = 0;
        int inverv;
        for (int i = 0; i < M; i++) {
            inverv = intervs(i);
            fT.segment(offset, inverv).setConstant(cT(i) / inverv);
            offset += inverv;
        }
        return;
    }

    static inline void mergeToCoarseGradT(const Eigen::VectorXi &intervs,
                                          Eigen::VectorXd &fineGdT) {
        int M = intervs.size();
        int offset = 0;
        int inverv;
        for (int i = 0; i < M; i++) {
            inverv = intervs(i);
            fineGdT(i) = fineGdT.segment(offset, inverv).mean();
            offset += inverv;
        }
        return;
    }

    static inline double objectiveFunc(void *ptrObj,     ///<[in] 输入参数结构体
                                       const double *x, ///<[in] 当前优化状态
                                       double *grad, ///<[in] 计算当前梯度值
                                       const int n  ///<[in] 优化变量数量 n
    ) {
        // 从GCOPTER类内获得一些变量
        QNMOPT &obj = *(QNMOPT *) ptrObj;
        // 初始化优化维度
        const int dimT = obj.dimFreeT;
        const int dimP = obj.dimFreeP;

        // 初始话优化参数
        const double rh = obj.rho;
        // 使用map将参数和double数组进行绑定
        Eigen::Map<const Eigen::VectorXd> t(x, dimT);
        Eigen::Map<const Eigen::VectorXd> p(x + dimT, dimP);
        // 使用map将梯度向量和输出double绑定
        Eigen::Map<Eigen::VectorXd> gradt(grad, dimT);
        // 时间梯度，尺寸等于多项式段数
        Eigen::VectorXd proxyGradT(obj.fineN);
        Eigen::Map<Eigen::VectorXd> gradp(grad + dimT, dimP);

        forwardT(t, obj.coarseT, obj.softT, obj.sumT, obj.c2dfm);
        splitToFineT(obj.coarseT, obj.intervals, obj.fineT);

        for (int i = 0; i < obj.freePts.cols(); i++) {
            obj.freePts.col(i).x() = p(3 * i);
            obj.freePts.col(i).y() = p(3 * i + 1);
            obj.freePts.col(i).z() = p(3 * i + 2);
            obj.innerP.col(i * 2) = obj.freePts.col(i);
        }

        double cost;
        // 生成轨迹并计算代价和梯度
        obj.minco.generate(obj.innerP, obj.fineT);
        obj.minco.evalTrajCostGrad(obj.cons,       /// cons = 8
                                   obj.vmax, obj.amax,  /// 最大速度的加速度
                                   obj.chi, cost,   ///位置约束的输出
                                   proxyGradT,/// 大致的时间分配和位置的梯度
                                   obj.gdInPs);
        cost += rh * obj.coarseT.sum();
        proxyGradT.array() += rh;
        mergeToCoarseGradT(obj.intervals, proxyGradT);
        addLayerTGrad(t, proxyGradT, obj.softT, obj.sumT, obj.c2dfm);
        int idx = 0;
        for (int i = 0; i < obj.gdInPs.cols(); i++) {
            if (i % 2 == 0) {
                Vec3 grd(obj.gdInPs.col(i).x(), obj.gdInPs.col(i).y(), obj.gdInPs.col(i).z());
                gradp.segment(idx * 3, 3) = grd;
                idx++;
            }
        }
        gradt = proxyGradT.head(dimT);
        return cost;

        //        addLayerPGrad(p, obj.idxVs, obj.cfgVs, gdFull, gradp, obj.retraction);
        //        fmt::print(fg(fmt::color::gold), " -- [DEBUG] cost = {} \n", cost);
        //        cout << "\t T vec:" << obj.fineT.transpose() << endl;
        //        cout << "\t gradT vec:" << proxyGradT.transpose() << endl;

        //        cout<<proxyGradT.transpose()<<endl;
        //        fmt::print(fg(fmt::color::gold)," -- [DEBUG] cost = {} \n",cost);
        //        mergeToCoarseGradT(obj.intervals, proxyGradT);
        //        addLayerTGrad(t, proxyGradT, obj.softT, obj.sumT, obj.c2dfm);
        //        addLayerPGrad(p, obj.idxVs, obj.cfgVs, obj.gdInPs, gradp, obj.retraction);
        //        if (obj.retraction) {
        //            restrictNorm(p, obj.idxVs, obj.cfgVs, cost, gradp);
        //        }
    }

public:
    typedef shared_ptr<QNMOPT> Ptr;

    /// gridMesh - 求栅格的分配
    inline void gridMesh(const Eigen::Matrix3d &iState, ///<[in] 输入起点状态
                         const Eigen::Matrix3d &fState,///<[in] 输入终点状态
                         const std::vector<Eigen::Matrix3Xd> &cfgPolyVs, ///<[in] 输入轨迹上的走廊顶点（每个走廊第一个点为坐标系）
                         const double &gridResolution,  ///<[out] 输出栅格分辨率
                         Eigen::VectorXi &intervalsVec) const///<[out] 输出时间分配
    {
        // M为走廊数量
        int M = intervalsVec.size();

        int curInterval, k;
        Eigen::Vector3d lastP, curP;
        // curP为起点位置
        curP = iState.col(0);
        // 迭代剩下M-1个点
        for (int i = 0; i < M - 1; i++) {
            // 更新上次位置
            lastP = curP;
            // k为第2*i+1个走廊的顶点坐标
            // 根据PolyVs的生成规律 [1 12 2 23 3 34]
            // 2 * i + 1都是交界处走廊的坐标
            k = cfgPolyVs[2 * i + 1].cols() - 1;
            // 将k个顶点求和，病除以（1+k） + 0 即可得到世界坐标系下的走廊中心
            curP = cfgPolyVs[2 * i + 1].rightCols(k).rowwise().sum() / (1.0 + k) +
                   cfgPolyVs[2 * i + 1].col(0);
            // 时间分配为INF, 所以这一段的时间分配为0，在加上ceil向上取整，这里就全部为1
            curInterval = ceil((curP - lastP).norm() / gridResolution);
            intervalsVec(i) = curInterval > 0 ? curInterval : 1;
        }
        lastP = curP;
        curP = fState.col(0);
        curInterval = ceil((curP - lastP).norm() / gridResolution);
        intervalsVec(M - 1) = curInterval > 0 ? curInterval : 1;

        return;
    }

    /*
     * 简要的函数说明文字
     *  @param [in] param1 参数1说明
     *  @param [out] param2 参数2说明
     *  @return 返回值说明
     */
    inline bool extractVs(const std::vector<Eigen::Matrix<double, 6, -1>> &hPs,
                          std::vector<Eigen::Matrix3Xd> &vPs) const {
        // M为走廊个数减一，可能是走廊交界的个数
        const int M = hPs.size() - 1;
        // vPs是一个2M+1的3x3矩阵 目前不知道是干啥的。
        vPs.clear();
        vPs.reserve(2 * M + 1);


        int nv;
        Eigen::Matrix<double, 6, -1> curIH;
        Eigen::Matrix3Xd curIV, curIOB;

        for (int i = 0; i < M; i++) {
            // 迭代i为走廊交界的个数, hPs为原始走廊，curIV为当前走廊的顶点 3xX
            if (!geoutils_snap::enumerateVs(hPs[i], curIV)) {
                return false;
            }
            // nv: number of Verties 顶点个数
            nv = curIV.cols();
            curIOB.resize(3, nv);
            // IOB的第一列等于第一个顶点
            curIOB.col(0) = curIV.col(0);
            // IOB的后续列等于第一个顶点坐标系下的作坐标
            curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
            // 将当前走廊放入vPs
            vPs.push_back(curIOB);
            // curIH 是6xX的，列数等于当前走廊的超平面个数+下一走廊的超平面个数，也就是把两个相邻走廊凑一起了
            curIH.resize(6, hPs[i].cols() + hPs[i + 1].cols());
            curIH.leftCols(hPs[i].cols()) = hPs[i];
            curIH.rightCols(hPs[i + 1].cols()) = hPs[i + 1];
            //继续找顶点
            if (!geoutils_snap::enumerateVs(curIH, curIV)) {
                return false;
            }
            nv = curIV.cols();
            curIOB.resize(3, nv);
            curIOB.col(0) = curIV.col(0);
            curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
            // 将相交走廊放进vPs
            vPs.push_back(curIOB);
        }
        // 至此，vPs中就存放了 [1 12 2 23 3 34 ... (M-1M)] 个走廊的顶点
        // 再放入最后一个走廊
        if (!geoutils_snap::enumerateVs(hPs.back(), curIV)) {
            return false;
        }
        nv = curIV.cols();
        curIOB.resize(3, nv);
        curIOB.col(0) = curIV.col(0);
        curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
        vPs.push_back(curIOB);

        return true;
    }

    inline bool setup(const double &rh,
                      const double &st,
                      const StatePVAJ &iniState,
                      const StatePVAJ &finState,
                      const vector<Vec3> waypts,
                      const double &gridRes,
                      const int &itgSpaces,
                      const double &vm,
                      const double &am,
                      const double &smoEps,
                      const Eigen::Vector3d &w,
                      bool c2diffeo,
                      bool retract) {
        // 是否使用c2连续可微
        c2dfm = c2diffeo;
        // 是否使用retraction，目前设置为否，汪博说是
        retraction = retract;
        // 是否使用soft时间，目前不知道什么意思，这只为true
        softT = rh > 0;
        if (softT) {
            rho = rh;
            sumT = 1.0;
        } else {
            rho = 0.0;
            sumT = st;
        }
        // 初始化边界条件
        iState = iniState;
        fState = finState;
        wayPts = waypts;
        // 需要优化的时间为waypts数量+1
        dimFreeT = (waypts.size() + 1) * 2;
        dimFreeP = (waypts.size() + 1) * 3;

        coarseN = (waypts.size() + 1) * 2;
        fineN = coarseN;
        smoothEps = smoEps;
        chi = w;
        vmax = vm;
        amax = am;
        cons.resize(fineN);
        cons.setConstant(itgSpaces);
        intervals.resize(dimFreeT);
        intervals.setOnes();
        // 将所有的初始状态限幅，避免无解
        double tempNorm;
        tempNorm = iState.col(1).norm();
        iState.col(1) *= tempNorm > vmax ? (vmax / tempNorm) : 1.0;
        tempNorm = fState.col(1).norm();
        fState.col(1) *= tempNorm > vmax ? (vmax / tempNorm) : 1.0;
        tempNorm = iState.col(2).norm();
        iState.col(2) *= tempNorm > amax ? (amax / tempNorm) : 1.0;
        tempNorm = fState.col(2).norm();
        fState.col(2) *= tempNorm > amax ? (amax / tempNorm) : 1.0;

        // Setup for L-BFGS solver
        // 初始话lbfgs的参数
        lbfgs::lbfgs_load_default_parameters(&lbfgs_params);

        // Allocate temp variables
        // 根据当前走廊的状态，初始化分配内存空间
        coarseT.resize(dimFreeT);
        fineT.resize(fineN);
        innerP.resize(3, fineN - 1);
        gdInPs.resize(3, fineN - 1);


        minco.reset(iniState, finState, smoothEps, fineN);

        return true;
    }

    inline void setInitial(const std::vector<Eigen::Matrix3Xd> &cfgPolyVs,///<[in] 飞行走廊顶点
                           const Eigen::VectorXi &intervs,  ///<[in] 有效段数
                           Eigen::VectorXd &vecT,   ///<[out] 初始时间分配
                           Eigen::Matrix3Xd &vecInP) const { ///<[out] 初始waypoint分配
        // 初始化时间分配直接拉满到速度
        const double allocationSpeed = vmax;
        int M = vecT.size();
        Eigen::Vector3d lastP, curP, delta;
        int offset, interv, k;

        offset = 0;
        curP = iState.col(0);
        for (int i = 0; i < M - 1; i++) {
            lastP = curP;
            interv = intervs(i);
            k = cfgPolyVs[2 * i + 1].cols() - 1;
            curP = cfgPolyVs[2 * i + 1].rightCols(k).rowwise().sum() / (1.0 + k) +
                   cfgPolyVs[2 * i + 1].col(0);
            delta = curP - lastP;
            vecT(i) = delta.norm() / allocationSpeed;
            delta /= interv;
            for (int j = 0; j < interv; j++) {
                vecInP.col(offset++) = (j + 1) * delta + lastP;
            }
        }
        interv = intervs(M - 1);
        lastP = curP;
        curP = fState.col(0);
        delta = curP - lastP;
        vecT(M - 1) = delta.norm() / allocationSpeed;
        delta /= interv;
        for (int j = 0; j < interv - 1; j++) {
            vecInP.col(offset++) = (j + 1) * delta + lastP;
        }

        return;
    }


    inline double optimize(Trajectory &traj,
                           double &relCostTol) {

/*        cout<<"iniState===============\n "<< iState<<endl;
        for(int i = 0 ; i < wayPts.size() ; i++){
            cout<< "Waypts "<<i<<" :"<<wayPts[i].transpose()<<endl;
        }g
        cout<<"finState===============\n "<< fState<<endl;*/

        //        fmt::print("OPT VAR NUM = {}\n",dimFreeT+dimFreeP);
        // 分配优化状态向量 尺寸等于优化时间和状态
        double *x = new double[dimFreeT + dimFreeP];
        // 使用map操作将时间优化与t绑定，状态优化与p绑定
        Eigen::Map<Eigen::VectorXd> t(x, dimFreeT);
        Eigen::Map<Eigen::VectorXd> p(x + dimFreeT, dimFreeP);
        // 初始化时间分配直接拉满到速度
        const double allocationSpeed = vmax;
        int M = dimFreeT;

        if (wayPts.size() == 0) {
            freePts.resize(3, 1);
            freePts = (fState.col(0) + iState.col(0)) / 2;
            innerP.col(0) = freePts;
            coarseT(0) = (fState.col(0) - iState.col(0)).norm() / allocationSpeed / 2;
            coarseT(1) = (fState.col(0) - iState.col(0)).norm() / allocationSpeed / 2;
        } else {
            freePts.resize(3, wayPts.size() + 1);
            freePts.col(0) = (wayPts[0] + iState.col(0)) / 2;
            for (int i = 0; i < wayPts.size() - 1; i++) {
                Vec3 midpt = (wayPts[i] + wayPts[i + 1]) / 2;
                freePts.col(i + 1) = midpt;
            }

            freePts.rightCols(1) = (wayPts[wayPts.size() - 1] + fState.col(0)) / 2;

            for (int i = 0; i < wayPts.size(); i++) {
                innerP.col(2 * (i)) = freePts.col(i);
                innerP.col(2 * (i) + 1) = wayPts[i];
            }

            innerP.rightCols(1) = freePts.rightCols(1);

            M = dimFreeT;
            Eigen::Vector3d lastP, curP, delta;
            curP = iState.col(0);
            for (int i = 0; i < M - 1; i++) {
                lastP = curP;
                curP = innerP.col(i);
                delta = curP - lastP;
                coarseT(i) = delta.norm() / allocationSpeed;
            }
            delta = innerP.rightCols(1) - fState.col(0);
            coarseT(dimFreeT - 1) = delta.norm() / allocationSpeed;
        }


        // 将优化器的时间初始化成一种类似log的形式
        backwardT(coarseT, t, softT, c2dfm);
        for (int i = 0; i < freePts.cols(); i++) {
            p(3 * i) = freePts.col(i).x();
            p(3 * i + 1) = freePts.col(i).y();
            p(3 * i + 2) = freePts.col(i).z();
        }

        // 初始主优化参数
        double minObjectivePenalty;
        lbfgs_params.mem_size = 64;
        lbfgs_params.past = 3;
        lbfgs_params.g_epsilon = 1.0e-16;
        lbfgs_params.min_step = 1.0e-32;
        lbfgs_params.abs_curv_cond = 0;
        lbfgs_params.delta = 1e-6;
        //        fmt::print(fg(fmt::color::gold), " -- [DEBUG] lbfgs init success.\n");
        //        fmt::print(fg(fmt::color::azure), "\tdimFreeT = {}\n", dimFreeT);
        // 开始下降
        int retCode = lbfgs::lbfgs_optimize(dimFreeT + dimFreeP,
                                            x,
                                            &minObjectivePenalty,
                                            &QNMOPT::objectiveFunc,
                                            nullptr,
                                            nullptr,
                                            this,
                                            &lbfgs_params);

        if (retCode < 0) {
            fmt::print(fg(fmt::color::red), " -- [DEBUG] lbfgs ret code = {}\n", retCode);
        }

        //        fineT = t.head(dimFreeT);
        //        for(int i = 0 ; i < dimFreeT ; i ++){
        //            fmt::print(" i-th fine T = {}\n", fineT(i));
        //        }
        forwardT(t, coarseT, softT, sumT, c2dfm);
        splitToFineT(coarseT, intervals, fineT);
        // 使用边界条件和时间向量生成一次minsnap

        minco.generate(innerP, fineT);
        traj = minco.getTraj();
        relCostTol = minObjectivePenalty;

        delete[] x;
        return minco.getTrajSnapCost();
    }
};

#endif