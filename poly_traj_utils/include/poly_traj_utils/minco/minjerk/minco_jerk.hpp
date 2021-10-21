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

#ifndef SE3GCOPTER_HPP
#define SE3GCOPTER_HPP


#include "poly_traj_utils/root_finder.hpp"
#include "lbfgs.hpp"
#include "./geoutils.hpp"
#include <Eigen/Eigen>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <vector>
#include <ctime>
#include "poly_traj_utils/poly_data_structure.hpp"
#include "poly_traj_utils/scope_timer.hpp"
#include <iomanip>      // std::setprecision

class MINCO_S3 {
public:
    MINCO_S3() = default;

    ~MINCO_S3() { A.destroy(); }
    //my
    // double compute_time = 0;

    inline void savelog2txt() {
        time_t t = time(NULL);
        char ch[64] = {0};
        char result[100] = {0};
        strftime(ch, sizeof(ch) - 1, "%Y%m%d--%H%M%S", localtime(&t));
        sprintf(result, "%s", ch);
        std::string filename = std::string(result);
        filename = "/home/yunfan/Logs/opt/opt.txt";// + filename + ".txt";
        std::ofstream log_file(filename);
        std::ostream_iterator<Eigen::Matrix<double, 1, -1>> log_iterator(log_file, "\n");
        std::copy(opt_log.begin(), opt_log.end(), log_iterator);
        printf(" -- [EXIT] \033[32m SAVE MPC LOGS SUCCESS! \033[0m\n");
    }

    inline void resetlog() {
        opt_log.clear();
    }

    inline void setSmootheps(double smo){
        smoothEps = smo;
    }

    inline void write2log(Eigen::VectorXd curlog) {
        opt_log.push_back(curlog.transpose());
    };
private:
    int N;
    Eigen::Matrix3d headPVA;
    Eigen::Matrix3d tailPVA;
    Eigen::VectorXd T1;
    BandedSystem A;
    Eigen::MatrixXd b;

    // Temp variables
    Eigen::VectorXd T2;
    Eigen::VectorXd T3;
    Eigen::VectorXd T4;
    Eigen::VectorXd T5;
    Eigen::MatrixXd gdC;
    double *cost_block;
    double smoothEps;

private:
    vector<Eigen::Matrix<double, 1, -1>> opt_log;

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


    template<typename EIGENVEC>
    inline void addGradJbyT(EIGENVEC &gdT) const {
        for (int i = 0; i < N; i++) {
            gdT(i) += 36.0 * b.row(6 * i + 3).squaredNorm() +
                      288.0 * b.row(6 * i + 4).dot(b.row(6 * i + 3)) * T1(i) +
                      576.0 * b.row(6 * i + 4).squaredNorm() * T2(i) +
                      720.0 * b.row(6 * i + 5).dot(b.row(6 * i + 3)) * T2(i) +
                      2880.0 * b.row(6 * i + 5).dot(b.row(6 * i + 4)) * T3(i) +
                      3600.0 * b.row(6 * i + 5).squaredNorm() * T4(i);
        }
        return;
    }

    template<typename EIGENMAT>
    inline void addGradJbyC(EIGENMAT &gdC) const {
        for (int i = 0; i < N; i++) {
            gdC.row(6 * i + 5) += 240.0 * b.row(6 * i + 3) * T3(i) +
                                  720.0 * b.row(6 * i + 4) * T4(i) +
                                  1440.0 * b.row(6 * i + 5) * T5(i);
            gdC.row(6 * i + 4) += 144.0 * b.row(6 * i + 3) * T2(i) +
                                  384.0 * b.row(6 * i + 4) * T3(i) +
                                  720.0 * b.row(6 * i + 5) * T4(i);
            gdC.row(6 * i + 3) += 72.0 * b.row(6 * i + 3) * T1(i) +
                                  144.0 * b.row(6 * i + 4) * T2(i) +
                                  240.0 * b.row(6 * i + 5) * T3(i);
        }
        return;
    }

    inline void solveAdjGradC(Eigen::MatrixXd &gdC) const {
        A.solveAdj(gdC);
        return;
    }

    template<typename EIGENVEC>
    inline void addPropCtoT(const Eigen::MatrixXd &adjGdC, EIGENVEC &gdT) const {
        Eigen::MatrixXd B1(6, 3), B2(3, 3);

        Eigen::RowVector3d negVel, negAcc, negJer, negSnp, negCrk;

        for (int i = 0; i < N - 1; i++) {
            negVel = -(b.row(i * 6 + 1) +
                       2.0 * T1(i) * b.row(i * 6 + 2) +
                       3.0 * T2(i) * b.row(i * 6 + 3) +
                       4.0 * T3(i) * b.row(i * 6 + 4) +
                       5.0 * T4(i) * b.row(i * 6 + 5));
            negAcc = -(2.0 * b.row(i * 6 + 2) +
                       6.0 * T1(i) * b.row(i * 6 + 3) +
                       12.0 * T2(i) * b.row(i * 6 + 4) +
                       20.0 * T3(i) * b.row(i * 6 + 5));
            negJer = -(6.0 * b.row(i * 6 + 3) +
                       24.0 * T1(i) * b.row(i * 6 + 4) +
                       60.0 * T2(i) * b.row(i * 6 + 5));
            negSnp = -(24.0 * b.row(i * 6 + 4) +
                       120.0 * T1(i) * b.row(i * 6 + 5));
            negCrk = -120.0 * b.row(i * 6 + 5);

            B1 << negSnp, negCrk, negVel, negVel, negAcc, negJer;

            gdT(i) += B1.cwiseProduct(adjGdC.block<6, 3>(6 * i + 3, 0)).sum();
        }

        negVel = -(b.row(6 * N - 5) +
                   2.0 * T1(N - 1) * b.row(6 * N - 4) +
                   3.0 * T2(N - 1) * b.row(6 * N - 3) +
                   4.0 * T3(N - 1) * b.row(6 * N - 2) +
                   5.0 * T4(N - 1) * b.row(6 * N - 1));
        negAcc = -(2.0 * b.row(6 * N - 4) +
                   6.0 * T1(N - 1) * b.row(6 * N - 3) +
                   12.0 * T2(N - 1) * b.row(6 * N - 2) +
                   20.0 * T3(N - 1) * b.row(6 * N - 1));
        negJer = -(6.0 * b.row(6 * N - 3) +
                   24.0 * T1(N - 1) * b.row(6 * N - 2) +
                   60.0 * T2(N - 1) * b.row(6 * N - 1));

        B2 << negVel, negAcc, negJer;

        gdT(N - 1) += B2.cwiseProduct(adjGdC.block<3, 3>(6 * N - 3, 0)).sum();

        return;
    }

    template<typename EIGENMAT>
    inline void addPropCtoP(const Eigen::MatrixXd &adjGdC, EIGENMAT &gdInP) const {
        for (int i = 0; i < N - 1; i++) {
            gdInP.col(i) += adjGdC.row(6 * i + 5).transpose();
        }
        return;
    }

    inline void normalizeFDF(const Eigen::Vector3d &x,
                             Eigen::Vector3d &xNor,
                             Eigen::Matrix3d &G) const {
        const double a = x(0), b = x(1), c = x(2);
        const double aSqr = a * a, bSqr = b * b, cSqr = c * c;
        const double ab = a * b, bc = b * c, ca = c * a;
        const double xSqrNorm = aSqr + bSqr + cSqr;
        const double xNorm = sqrt(xSqrNorm);
        const double den = xSqrNorm * xNorm;
        xNor = x / xNorm;
        G(0, 0) = bSqr + cSqr;
        G(0, 1) = -ab;
        G(0, 2) = -ca;
        G(1, 0) = -ab;
        G(1, 1) = aSqr + cSqr;
        G(1, 2) = -bc;
        G(2, 0) = -ca;
        G(2, 1) = -bc;
        G(2, 2) = aSqr + bSqr;
        G /= den;
        return;
    }

    template<typename EIGENVEC>
    inline void addTimeIntPenalty(const Eigen::VectorXi cons,
                                  const double vMax,
                                  const double aMax,
                                  const double thrAccMin,
                                  const double thrAccMax,
                                  const double bdrMax,
                                  const double gAcc,
                                  const Eigen::Vector4d ci,
                                  double &cost,
                                  EIGENVEC &gdT,
                                  Eigen::MatrixXd &gdC) const {
        double pena = 0.0;
        const double vMaxSqr = vMax * vMax;
        const double aMaxSqr = aMax * aMax;
        const double thrAccMinSqr = thrAccMin * thrAccMin;
        const double thrAccMaxSqr = thrAccMax * thrAccMax;
        const double bdrMaxSqr = bdrMax * bdrMax;

        Eigen::Vector3d pos, vel, acc, jer, sna;
        double step, alpha;
        double s1, s2, s3, s4, s5;
        Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3, beta4;

        Eigen::Matrix3d rotM;

        Eigen::Vector3d h, zB, czB, xC, yB, xB;
        Eigen::Matrix3d dnczB, dzB, cdzB, dyB, dxB;
        Eigen::Vector3d outerNormal, point;

        Eigen::Vector3d eNormGd;
        Eigen::Matrix3d gradSdTxyz;
        Eigen::Vector3d gradSdT;

        Eigen::Matrix<double, 6, 3> gradSdCx, gradSdCy, gradSdCz, gradSdC;
        Eigen::Matrix<double, 6, 3> beta2dOuterNormalTp, beta0dOuterNormalTp;

        double violaVel, violaAcc, violaThrl, violaThrh, violaBdr;
        double violaVelPenaD, violaAccPenaD, violaThrlPenaD, violaThrhPenaD, violaBdrPenaD;
        double violaVelPena, violaAccPena, violaThrlPena, violaThrhPena, violaBdrPena;
        Eigen::Matrix<double, 6, 3> gradViolaVc, gradViolaAc, gradViolaThrlc, gradViolaThrhc, gradViolaBdrc;
        double gradViolaVt,gradViolaAt, gradViolaThrlt, gradViolaThrht, gradViolaBdrt;
        double fThr, sqrMagThr, sqrMagBdr;
        Eigen::Vector3d dfThr, dSqrMagThr, bdr, xyBdr;
        Eigen::Vector3d dSqrMagBdr, rotTrDotJer;
        Eigen::Matrix3d dBdr, dxyBdr;
        Eigen::Vector3d dJerSqrMagBdr;
        Eigen::Matrix3d dJerBdr, dJerxyBdr;
        double omg;

        int innerLoop;
        for (int i = 0; i < N; i++) {
            const auto &c = b.block<6, 3>(i * 6, 0);
            step = T1(i) / cons(i);
            s1 = 0.0;
            innerLoop = cons(i) + 1;
            for (int j = 0; j < innerLoop; j++) {
                s2 = s1 * s1;
                s3 = s2 * s1;
                s4 = s2 * s2;
                s5 = s4 * s1;
                beta0 << 1.0, s1, s2, s3, s4, s5;
                beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4;
                beta2 << 0.0, 0.0, 2.0, 6.0 * s1, 12.0 * s2, 20.0 * s3;
                beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * s1, 60.0 * s2;
                beta4 << 0.0, 0.0, 0.0, 0.0, 24.0, 120.0 * s1;
                alpha = 1.0 / cons(i) * j;

                vel = c.transpose() * beta1;
                acc = c.transpose() * beta2;
                jer = c.transpose() * beta3;
                sna = c.transpose() * beta4;

                if(bdrMax<0 && thrAccMax < 0 && thrAccMin < 0 ){

                    violaVel = vel.squaredNorm() - vMaxSqr;
                    violaAcc = acc.squaredNorm() - aMaxSqr;
                    omg = (j == 0 || j == innerLoop - 1) ? 0.5 : 1.0;

                    if (violaVel > 0.0) {
                        positiveSmoothedL1(violaVel, violaVelPena, violaVelPenaD);
                        gradViolaVc = 2.0 * beta1 * vel.transpose();
                        gradViolaVt = 2.0 * alpha * vel.transpose() * acc;
                        gdC.block<6, 3>(i * 6, 0) += omg * step * ci(1) * violaVelPenaD * gradViolaVc;
                        gdT(i) += omg * (ci(1) * violaVelPenaD * gradViolaVt * step +
                                ci(1) * violaVelPena / cons(i));
                        pena += omg * step * ci(1) * violaVelPena;
                    }

                    if (violaAcc > 0.0) {
                        positiveSmoothedL1(violaAcc, violaAccPena, violaAccPenaD);
                        gradViolaAc = 2.0 * beta2 * acc.transpose();
                        gradViolaAt = 2.0 * alpha * acc.transpose() * jer;
                        gdC.block<6, 3>(i * 6, 0) += omg * step * ci(2) * violaAccPenaD * gradViolaAc;
                        gdT(i) += omg * (ci(2) * violaAccPenaD * gradViolaAt * step +
                                ci(2) * violaAccPena / cons(i));
                        pena += omg * step * ci(2) * violaAccPena;
                    }

                }else{

                    h = acc;
                    h(2) += gAcc;
                    normalizeFDF(h, zB, dzB);
                    //Zb
                    czB << 0.0, zB(2), -zB(1);
                    cdzB << Eigen::RowVector3d::Zero(), dzB.row(2), -dzB.row(1);
                    normalizeFDF(czB, yB, dnczB);
                    //Yb
                    xB = yB.cross(zB);
                    dyB = dnczB * cdzB;
                    dxB.col(0) = dyB.col(0).cross(zB) + yB.cross(dzB.col(0));
                    dxB.col(1) = dyB.col(1).cross(zB) + yB.cross(dzB.col(1));
                    dxB.col(2) = dyB.col(2).cross(zB) + yB.cross(dzB.col(2));
                    rotM << xB, yB, zB;

                    gradSdTxyz.col(0) = dxB * jer;
                    gradSdTxyz.col(1) = dyB * jer;
                    gradSdTxyz.col(2) = dzB * jer;

                    fThr = h.norm();
                    dfThr = h / fThr;
                    sqrMagThr = fThr * fThr;
                    dSqrMagThr = 2 * h;
                    rotTrDotJer = rotM.transpose() * jer;
                    bdr = rotTrDotJer / fThr;
                    xyBdr << -bdr(1), bdr(0), 0.0;
                    sqrMagBdr = xyBdr.squaredNorm();
                    dBdr = -rotTrDotJer * dfThr.transpose() / (fThr * fThr) -
                            rotM.transpose() * (dxB * rotTrDotJer(0) + dyB * rotTrDotJer(1) + dzB * rotTrDotJer(2)) / fThr;
                    dxyBdr << -dBdr.row(1), dBdr.row(0), Eigen::RowVector3d::Zero();
                    dSqrMagBdr = 2.0 * xyBdr.transpose() * dxyBdr;
                    dJerBdr = rotM.transpose() / fThr;
                    dJerxyBdr << -dJerBdr.row(1), dJerBdr.row(0), Eigen::RowVector3d::Zero();
                    dJerSqrMagBdr = 2.0 * xyBdr.transpose() * dJerxyBdr;

                    violaVel = vel.squaredNorm() - vMaxSqr;
                    violaThrl = thrAccMinSqr - sqrMagThr;
                    violaThrh = sqrMagThr - thrAccMaxSqr;
                    violaBdr = sqrMagBdr - bdrMaxSqr;

                    omg = (j == 0 || j == innerLoop - 1) ? 0.5 : 1.0;


                    if (violaThrl > 0.0) {
                        positiveSmoothedL1(violaThrl, violaThrlPena, violaThrlPenaD);
                        gradViolaThrlc = -beta2 * dSqrMagThr.transpose();
                        gradViolaThrlt = -alpha * dSqrMagThr.transpose() * jer;
                        gdC.block<6, 3>(i * 6, 0) += omg * step * ci(2) * violaThrlPenaD * gradViolaThrlc;
                        gdT(i) += omg * (ci(2) * violaThrlPenaD * gradViolaThrlt * step +
                                ci(2) * violaThrlPena / cons(i));
                        pena += omg * step * ci(2) * violaThrlPena;
                    }

                    if (violaThrh > 0.0) {
                        positiveSmoothedL1(violaThrh, violaThrhPena, violaThrhPenaD);
                        gradViolaThrhc = beta2 * dSqrMagThr.transpose();
                        gradViolaThrht = alpha * dSqrMagThr.transpose() * jer;
                        gdC.block<6, 3>(i * 6, 0) += omg * step * ci(2) * violaThrhPenaD * gradViolaThrhc;
                        gdT(i) += omg * (ci(2) * violaThrhPenaD * gradViolaThrht * step +
                                ci(2) * violaThrhPena / cons(i));
                        pena += omg * step * ci(2) * violaThrhPena;
                    }

                    if (violaBdr > 0.0) {
                        positiveSmoothedL1(violaBdr, violaBdrPena, violaBdrPenaD);
                        gradViolaBdrc = beta2 * dSqrMagBdr.transpose() + beta3 * dJerSqrMagBdr.transpose();
                        gradViolaBdrt = alpha * (dSqrMagBdr.dot(jer) + dJerSqrMagBdr.dot(sna));
                        gdC.block<6, 3>(i * 6, 0) += omg * step * ci(3) * violaBdrPenaD * gradViolaBdrc;
                        gdT(i) += omg * (ci(3) * violaBdrPenaD * gradViolaBdrt * step +
                                ci(3) * violaBdrPena / cons(i));
                        pena += omg * step * ci(3) * violaBdrPena;
                    }

                    if (violaVel > 0.0) {
                        positiveSmoothedL1(violaVel, violaVelPena, violaVelPenaD);
                        gradViolaVc = 2.0 * beta1 * vel.transpose();
                        gradViolaVt = 2.0 * alpha * vel.transpose() * acc;
                        gdC.block<6, 3>(i * 6, 0) += omg * step * ci(1) * violaVelPenaD * gradViolaVc;
                        gdT(i) += omg * (ci(1) * violaVelPenaD * gradViolaVt * step +
                                ci(1) * violaVelPena / cons(i));
                        pena += omg * step * ci(1) * violaVelPena;

                    }

                }
                s1 += step;
            }
        }
        // end_time = clock();
        // double endtime=(double)(end_time-start_time)/CLOCKS_PER_SEC;
        // compute_time+=endtime;
        cost += pena;
        return;
    }

public:
    inline void reset(const Eigen::Matrix3d &headState,
                      const Eigen::Matrix3d &tailState,
                      const int &pieceNum) {
        N = pieceNum;
        headPVA = headState;
        tailPVA = tailState;
        T1.resize(N);
        A.create(6 * N, 6, 6);
        b.resize(6 * N, 3);
        gdC.resize(6 * N, 3);
        return;
    }

    inline void generate(const Eigen::MatrixXd &inPs,
                         const Eigen::VectorXd &ts) {
        T1 = ts;
        T2 = T1.cwiseProduct(T1);
        T3 = T2.cwiseProduct(T1);
        T4 = T2.cwiseProduct(T2);
        T5 = T4.cwiseProduct(T1);

        A.reset();
        b.setZero();

        A(0, 0) = 1.0;
        A(1, 1) = 1.0;
        A(2, 2) = 2.0;
        b.row(0) = headPVA.col(0).transpose();
        b.row(1) = headPVA.col(1).transpose();
        b.row(2) = headPVA.col(2).transpose();

        for (int i = 0; i < N - 1; i++) {
            A(6 * i + 3, 6 * i + 3) = 6.0;
            A(6 * i + 3, 6 * i + 4) = 24.0 * T1(i);
            A(6 * i + 3, 6 * i + 5) = 60.0 * T2(i);
            A(6 * i + 3, 6 * i + 9) = -6.0;
            A(6 * i + 4, 6 * i + 4) = 24.0;
            A(6 * i + 4, 6 * i + 5) = 120.0 * T1(i);
            A(6 * i + 4, 6 * i + 10) = -24.0;
            A(6 * i + 5, 6 * i) = 1.0;
            A(6 * i + 5, 6 * i + 1) = T1(i);
            A(6 * i + 5, 6 * i + 2) = T2(i);
            A(6 * i + 5, 6 * i + 3) = T3(i);
            A(6 * i + 5, 6 * i + 4) = T4(i);
            A(6 * i + 5, 6 * i + 5) = T5(i);
            A(6 * i + 6, 6 * i) = 1.0;
            A(6 * i + 6, 6 * i + 1) = T1(i);
            A(6 * i + 6, 6 * i + 2) = T2(i);
            A(6 * i + 6, 6 * i + 3) = T3(i);
            A(6 * i + 6, 6 * i + 4) = T4(i);
            A(6 * i + 6, 6 * i + 5) = T5(i);
            A(6 * i + 6, 6 * i + 6) = -1.0;
            A(6 * i + 7, 6 * i + 1) = 1.0;
            A(6 * i + 7, 6 * i + 2) = 2 * T1(i);
            A(6 * i + 7, 6 * i + 3) = 3 * T2(i);
            A(6 * i + 7, 6 * i + 4) = 4 * T3(i);
            A(6 * i + 7, 6 * i + 5) = 5 * T4(i);
            A(6 * i + 7, 6 * i + 7) = -1.0;
            A(6 * i + 8, 6 * i + 2) = 2.0;
            A(6 * i + 8, 6 * i + 3) = 6 * T1(i);
            A(6 * i + 8, 6 * i + 4) = 12 * T2(i);
            A(6 * i + 8, 6 * i + 5) = 20 * T3(i);
            A(6 * i + 8, 6 * i + 8) = -2.0;

            b.row(6 * i + 5) = inPs.col(i).transpose();
        }

        A(6 * N - 3, 6 * N - 6) = 1.0;
        A(6 * N - 3, 6 * N - 5) = T1(N - 1);
        A(6 * N - 3, 6 * N - 4) = T2(N - 1);
        A(6 * N - 3, 6 * N - 3) = T3(N - 1);
        A(6 * N - 3, 6 * N - 2) = T4(N - 1);
        A(6 * N - 3, 6 * N - 1) = T5(N - 1);
        A(6 * N - 2, 6 * N - 5) = 1.0;
        A(6 * N - 2, 6 * N - 4) = 2 * T1(N - 1);
        A(6 * N - 2, 6 * N - 3) = 3 * T2(N - 1);
        A(6 * N - 2, 6 * N - 2) = 4 * T3(N - 1);
        A(6 * N - 2, 6 * N - 1) = 5 * T4(N - 1);
        A(6 * N - 1, 6 * N - 4) = 2;
        A(6 * N - 1, 6 * N - 3) = 6 * T1(N - 1);
        A(6 * N - 1, 6 * N - 2) = 12 * T2(N - 1);
        A(6 * N - 1, 6 * N - 1) = 20 * T3(N - 1);

        b.row(6 * N - 3) = tailPVA.col(0).transpose();
        b.row(6 * N - 2) = tailPVA.col(1).transpose();
        b.row(6 * N - 1) = tailPVA.col(2).transpose();

        A.factorizeLU();
        A.solve(b);

        return;
    }

    inline double getTrajJerkCost() const {
        double objective = 0.0;
        for (int i = 0; i < N; i++) {
            objective += 36.0 * b.row(6 * i + 3).squaredNorm() * T1(i) +
                         144.0 * b.row(6 * i + 4).dot(b.row(6 * i + 3)) * T2(i) +
                         192.0 * b.row(6 * i + 4).squaredNorm() * T3(i) +
                         240.0 * b.row(6 * i + 5).dot(b.row(6 * i + 3)) * T3(i) +
                         720.0 * b.row(6 * i + 5).dot(b.row(6 * i + 4)) * T4(i) +
                         720.0 * b.row(6 * i + 5).squaredNorm() * T5(i);
        }
        return objective;
    }

    template<typename EIGENVEC, typename EIGENMAT>
    inline void evalTrajCostGrad(const Eigen::VectorXi &cons,
                                 const double &vMax,
                                 const double &aMax,
                                 const double &thrAccMin,
                                 const double &thrAccMax,
                                 const double &bdrMax,
                                 const double &gAcc,
                                 const Eigen::Vector4d &ci,
                                 double &cost,
                                 EIGENVEC &gdT,
                                 EIGENMAT &gdInPs) {
//        TimeConsuming t__("evalTrajCostGrad");
        gdT.setZero();
        gdInPs.setZero();
        gdC.setZero();

        cost = getTrajJerkCost();
        addGradJbyT(gdT);
        addGradJbyC(gdC);

        addTimeIntPenalty(cons, vMax, aMax, thrAccMin, thrAccMax, bdrMax,
                          gAcc, ci, cost, gdT, gdC);

        solveAdjGradC(gdC);
        addPropCtoT(gdC, gdT);
        addPropCtoP(gdC, gdInPs);
    }

    inline Trajectory getTraj(void) const {
        Trajectory traj;
        traj.reserve(N);
        for (int i = 0; i < N; i++) {
            traj.emplace_back(T1(i), b.block<6, 3>(6 * i, 0).transpose().rowwise().reverse());
        }
        return traj;
    }
};

class SE3GCOPTER {
public:
    int cnt = 0;
    inline void cntpp(){
        cnt++;
    }

private:
    // Use C2 or Cinf diffeo
    bool c2dfm;

    // Use soft time or not
    bool softT;

    // Weight for time regularization term
    double rho;

    // Fixed total time
    double sumT;

    //Minimum Jerk Optimizer
    MINCO_S3 jerkOpt;

    // Temp variables for problem solving
    StatePVA iState;
    StatePVA fState;
    vector<Vec3> wayPts;


    // Each col of cfgHs denotes a facet (outter_normal^T,point^T)^T
    std::vector<Eigen::MatrixXd> cfgVs;
    std::vector<Eigen::MatrixXd> cfgHs;
    Eigen::MatrixXd gdInPs;

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
    Eigen::VectorXd freeT;
    Eigen::VectorXd coarseT;
    Eigen::VectorXd fineT;
    Eigen::MatrixXd innerP;
    Eigen::MatrixXd init_innerP;
    Eigen::Matrix3Xd freePts;
    // Params for constraints
    Eigen::VectorXi cons;
    Eigen::Vector4d chi;

    Eigen::Vector3d ellipsoid;
    double safeMargin;
    double vMax, aMax;
    double thrAccMin;
    double thrAccMax;
    double bdrMax;
    double gAcc;

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
    static inline void backwardT(const Eigen::VectorXd &vecT,
                                 EIGENVEC &t,
                                 bool soft,
                                 bool c2) {
        if (soft) {
            if (c2) {
                int M = vecT.size();
                for (int i = 0; i < M; i++) {
                    t(i) = vecT(i) > 1.0
                           ? (sqrt(2.0 * vecT(i) - 1.0) - 1.0)
                           : (1.0 - sqrt(2.0 / vecT(i) - 1.0));
                    /*
                    
                    t  = 2/[(tau-1)^2+1] tau<0
                       = 1/2[(tau+1)^2+1] tau>=0
                    */
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

    template<typename EIGENVEC>
    static inline void forwardP(const EIGENVEC &p,
                                const Eigen::VectorXi &idVs,
                                const std::vector<Eigen::MatrixXd> &cfgPolyVs,
                                Eigen::MatrixXd &inP) {
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
        return;
    }

    static inline double objectiveNLS(void *ptrPOBs,
                                      const double *x,
                                      double *grad,
                                      const int n) {
        const Eigen::MatrixXd &pobs = *(Eigen::MatrixXd *) ptrPOBs;
        Eigen::Map<const Eigen::VectorXd> p(x, n);
        Eigen::Map<Eigen::VectorXd> gradp(grad, n);

        double qnsqr = p.squaredNorm();
        double qnsqrp1 = qnsqr + 1.0;
        double qnsqrp1sqr = qnsqrp1 * qnsqrp1;
        Eigen::VectorXd r = 2.0 / qnsqrp1 * p;

        Eigen::Vector3d delta = pobs.rightCols(n) * r.cwiseProduct(r) +
                                pobs.col(1) - pobs.col(0);
        double cost = delta.squaredNorm();
        Eigen::Vector3d gradR3 = 2 * delta;

        Eigen::VectorXd gdr = pobs.rightCols(n).transpose() * gradR3;
        gdr = gdr.array() * r.array() * 2.0;
        gradp = gdr * 2.0 / qnsqrp1 -
                p * 4.0 * gdr.dot(p) / qnsqrp1sqr;

        return cost;
    }

    template<typename EIGENVEC>
    static inline void backwardP(const Eigen::MatrixXd &inP,
                                 const Eigen::VectorXi &idVs,
                                 const std::vector<Eigen::MatrixXd> &cfgPolyVs,
                                 EIGENVEC &p) {
        int M = inP.cols();
        int j = 0, k, idx;

        // Parameters for tiny nonlinear least squares
        double minSqrD;
        lbfgs::lbfgs_parameter_t nls_params;
        lbfgs::lbfgs_load_default_parameters(&nls_params);
        nls_params.g_epsilon = FLT_EPSILON;
        nls_params.max_iterations = 128;

        Eigen::MatrixXd pobs;
        for (int i = 0; i < M; i++) {
            idx = idVs(i);
            k = cfgPolyVs[idx].cols() - 1;
            p.segment(j, k).setConstant(1.0 / (sqrt(k + 1.0) + 1.0));
            pobs.resize(3, k + 2);
            pobs << inP.col(i), cfgPolyVs[idx];
            lbfgs::lbfgs_optimize(k,
                                  p.data() + j,
                                  &minSqrD,
                                  &SE3GCOPTER::objectiveNLS,
                                  nullptr,
                                  nullptr,
                                  &pobs,
                                  &nls_params);

            j += k;
        }

        return;
    }

    template<typename EIGENVEC>
    static inline void addLayerTGrad(const Eigen::VectorXd &t,
                                     EIGENVEC &gradT,
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

    template<typename EIGENVEC_0, typename EIGENVEC_1>
    static inline void addLayerPGrad(EIGENVEC_0 &p,
                                     const Eigen::VectorXi &idVs,
                                     const std::vector<Eigen::MatrixXd> &cfgPolyVs,
                                     const Eigen::MatrixXd &gradInPs,
                                     EIGENVEC_1 &grad) {
        int M = gradInPs.cols();

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


    static inline double objectiveFunc(void *ptrObj,
                                       const double *x,
                                       double *grad,
                                       const int n) {
//        TimeConsuming t__("objectiveFunc");
        SE3GCOPTER &obj = *(SE3GCOPTER *) ptrObj;
        const int dimT = obj.dimFreeT;
        const int dimP = obj.dimFreeP;
        const double rh = obj.rho;
        Eigen::Map<const Eigen::VectorXd> t(x, dimT);
        Eigen::Map<const Eigen::VectorXd> p(x + dimT, dimP);
        Eigen::Map<Eigen::VectorXd> gradt(grad, dimT);
        Eigen::VectorXd proxyGradT(obj.dimFreeT);
        Eigen::Map<Eigen::VectorXd> gradp(grad + dimT, dimP);

        // 从微分同胚的tau映射为T
        forwardT(t, obj.freeT, obj.softT, obj.sumT, obj.c2dfm);
//        splitToFineT(obj.coarseT, obj.intervals, obj.fineT);
//        forwardP(p, obj.idxVs, obj.cfgVs, obj.innerP);
        // 从优化变量中取出新的innerP
        for (int i = 0; i < obj.freePts.cols(); i++) {
            obj.freePts.col(i).x() = p(3 * i);
            obj.freePts.col(i).y() = p(3 * i + 1);
            obj.freePts.col(i).z() = p(3 * i + 2);
            obj.innerP.col(i * 2) = obj.freePts.col(i);
        }

        double cost;
        obj.jerkOpt.generate(obj.innerP, obj.freeT);
        obj.jerkOpt.evalTrajCostGrad(obj.cons, obj.vMax, obj.aMax,obj.thrAccMin,
                                     obj.thrAccMax, obj.bdrMax, obj.gAcc, obj.chi,
                                     cost, proxyGradT, obj.gdInPs);

        cost += rh * obj.freeT.sum();
        proxyGradT.array() += rh;

        addLayerTGrad(t, proxyGradT, obj.softT, obj.sumT, obj.c2dfm);

        // 将所有内部点的梯度信息，提取出freep的梯度加入到gradp中
        int idx = 0;
        for (int i = 0; i < obj.gdInPs.cols(); i++) {
            if (i % 2 == 0) {
                Vec3 grd(obj.gdInPs.col(i).x(), obj.gdInPs.col(i).y(), obj.gdInPs.col(i).z());
                gradp.segment(idx * 3, 3) = grd;
                idx++;
            }
        }
        gradt = proxyGradT.head(dimT);

//        cout<<std::setprecision(5)<< "it = " << iter_num++<<" gradT = "<<gradt.transpose()<<"\tgradP = "<<gradp.transpose()<<" freep = "<<p.transpose()<<" freet = "<<t.transpose()<<endl;

        obj.cntpp();
        return cost;
    }

public:

    inline bool setup(const double &rh,
                      const double &st,
                      const double &smo,
                      const StatePVA &iniState,
                      const StatePVA &finState,
                      const vector<Vec3> waypts,
                      const int &itgSpaces,
                      const double &vm,
                      const double &am,
                      const double &minThrAcc,
                      const double &maxThrAcc,
                      const double &bodyRateMax,
                      const double &g,
                      const Eigen::Vector4d &w,
                      bool c2diffeo) {
        // 是否使用c2连续可微
        c2dfm = c2diffeo;
        // 是否使用soft时间，目前不知道什么意思，这只为true
        softT = rh > 0;
        if (softT) {
            rho = rh;
            sumT = 1.0;
        } else {
            rho = 0.0;
            sumT = st;
        }
        // 初始化固定约束条件
        iState = iniState;
        fState = finState;
        wayPts = waypts;
        jerkOpt.setSmootheps(smo);
        // 需要优化的时间为waypts数量+1
        dimFreeT = (waypts.size() + 1) * 2;
        dimFreeP = (waypts.size() + 1) * 3;
        coarseN = (waypts.size() + 1) * 2;
        freeT.resize(dimFreeT);

        chi = w;
        vMax = vm;
        aMax = am;
        thrAccMin = minThrAcc;
        thrAccMax = maxThrAcc;
        bdrMax = bodyRateMax;
        gAcc = g;

        cons.resize(dimFreeT);
        cons.setConstant(itgSpaces);

        // 将初始速度进行保方向饱和
        double tempNorm;
        tempNorm = iState.col(1).norm();
        iState.col(1) *= tempNorm > vMax ? (vMax / tempNorm) : 1.0;
        tempNorm = fState.col(1).norm();
        fState.col(1) *= tempNorm > vMax ? (vMax / tempNorm) : 1.0;


        // Setup for L-BFGS solver
        // 初始化lbfgs的参数
        lbfgs::lbfgs_load_default_parameters(&lbfgs_params);


        // 初始化内部waypoint矩阵尺寸
        innerP.resize(3, 2 * waypts.size() + 1);
        gdInPs.resize(3, 2 * waypts.size() + 1);

        // 初始话Minsnap线性方程组尺寸
        jerkOpt.reset(iniState, finState, dimFreeT);

        return true;
    }

    inline void setInitial() {
        // 初始化时间分配直接拉满到速度
        const double allocationSpeed = vMax;

        if (wayPts.size() == 0) {
            // 如果只固定起点和终点，则只有一个自由waypoint
            freePts.resize(3, 1);
            // 将自由点设置为起点和终点的中间
            freePts = (fState.col(0) + iState.col(0)) / 2;
            // 轨迹内部所有waypoint等于自由点
            innerP.col(0) = freePts;
            // 时间分配为最大速度跑直线
            freeT(0) = (fState.col(0) - iState.col(0)).norm() / 2 / allocationSpeed;
            freeT(1) = (fState.col(0) - iState.col(0)).norm() / 2 / allocationSpeed;
        } else {
            // 如果有固定的waypoint了
            freePts.resize(3, wayPts.size() + 1);
            // 第一个自由点的坐标为第一个waypoint和起点中间
            freePts.col(0) = (wayPts[0] + iState.col(0)) / 2;
            // 中间的自由点坐标为两个相邻waypoint。
            for (int i = 0; i < wayPts.size() - 1; i++) {
                Vec3 midpt = (wayPts[i] + wayPts[i + 1]) / 2;
                freePts.col(i + 1) = midpt;
            }

            freePts.rightCols(1) = (wayPts[wayPts.size() - 1] + fState.col(0)) / 2;

            // 随后填充所有内部点。
            for (int i = 0; i < wayPts.size(); i++) {
                innerP.col(2 * (i)) = freePts.col(i);
                innerP.col(2 * (i) + 1) = wayPts[i];
            }

            innerP.rightCols(1) = freePts.rightCols(1);

            // 计算时间分配
            Eigen::Vector3d lastP, curP, delta;
            curP = iState.col(0);
            for (int i = 0; i < dimFreeT - 1; i++) {
                lastP = curP;
                curP = innerP.col(i);
                delta = curP - lastP;
                freeT(i) = delta.norm() / allocationSpeed;
            }
            delta = innerP.rightCols(1) - fState.col(0);
            freeT(dimFreeT - 1) = delta.norm() / allocationSpeed;
        }
        init_innerP = innerP;

        return;
    }

    inline vector<Vec3> getInnerPts() {
        vector<Vec3> inp;
        for (int i = 0; i < innerP.cols(); i++) {
            inp.push_back(innerP.col(i));
        }
        return inp;
    }

    inline vector<Vec3> getInitInnerPts() {
        vector<Vec3> inp;
        for (int i = 0; i < init_innerP.cols(); i++) {
            inp.push_back(init_innerP.col(i));
        }
        return inp;
    }

    inline double optimize(Trajectory &traj,
                           double &relCostTol) {
        TimeConsuming t__("Total opt time");
        // 分配优化状态向量 尺寸等于优化时间和状态
        double *x = new double[dimFreeT + dimFreeP];
        // 使用map操作将时间优化与t绑定，状态优化与p绑定
        Eigen::Map<Eigen::VectorXd> t(x, dimFreeT);
        Eigen::Map<Eigen::VectorXd> p(x + dimFreeT, dimFreeP);
        cnt = 0;
        // 初始化优化状态变量
        setInitial();

        // 将优化器时间填充为微分同胚的形式
        backwardT(freeT, t, softT, c2dfm);

        // 将优化状态直接填入优化变量
        for (int i = 0; i < freePts.cols(); i++) {
            p(3 * i) = freePts.col(i).x();
            p(3 * i + 1) = freePts.col(i).y();
            p(3 * i + 2) = freePts.col(i).z();
        }
        jerkOpt.resetlog();
        // 初始主优化参数
        double minObjectivePenalty;
        lbfgs_params.mem_size = 64;
        lbfgs_params.past = 3;
        lbfgs_params.g_epsilon = 1.0e-16;
        lbfgs_params.min_step = 1.0e-32;
        lbfgs_params.abs_curv_cond = 1;
        lbfgs_params.delta = 1e-4;
        //        fmt::print(fg(fmt::color::gold), " -- [DEBUG] lbfgs init success.\n");
        //        fmt::print(fg(fmt::color::azure), "\tdimFreeT = {}\n", dimFreeT);
        // 开始下降
        int retCode = lbfgs::lbfgs_optimize(dimFreeT + dimFreeP,
                                            x,
                                            &minObjectivePenalty,
                                            &SE3GCOPTER::objectiveFunc,
                                            nullptr,
                                            nullptr,
                                            this,
                                            &lbfgs_params);

        if (retCode < 0) {
            fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, " -- RETCODE {}. Backend failed, return!\n", retCode);
            cout << "iniState===============\n " << iState << endl;
            cout << "=================================" << endl;
            for (int i = 0; i < wayPts.size(); i++) {
                cout << "Waypts " << i << " :" << wayPts[i].transpose() << endl;
            }
            cout << "finState===============\n " << fState << endl;

        }
        else{
            fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, " -- RETCODE {}. Success with it_num = {}\n", retCode,cnt);
        }

        //        fineT = t.head(dimFreeT);
        //        for(int i = 0 ; i < dimFreeT ; i ++){
        //            fmt::print(" i-th fine T = {}\n", fineT(i));
        //        }
        // 从最终优化的时间中提取出真实的时间分配
        forwardT(t, freeT, softT, sumT, c2dfm);

        // 使用边界条件和时间向量生成一次minsnap
        jerkOpt.generate(innerP, freeT);
        traj = jerkOpt.getTraj();
        relCostTol = minObjectivePenalty;
//        jerkOpt.savelog2txt();
        fmt::print("dimT = {}, dimP ={}.\n", dimFreeT, dimFreeP);
        delete[] x;
        return jerkOpt.getTrajJerkCost();
    }
};

#endif
