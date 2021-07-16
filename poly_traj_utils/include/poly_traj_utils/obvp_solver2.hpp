//
// Created by yunfan on 2021/3/19.
//

#ifndef SRC_OBVP_SOLVER_HPP
#define SRC_OBVP_SOLVER_HPP

#include "root_finder.hpp"
#include "traj_utils.hpp"
#include "Eigen/Dense"
#include "vector"
#include "poly_visual_utils.hpp"

#include "scope_timer.hpp"

using namespace Eigen;
using namespace std;

typedef Eigen::Matrix<double, 12, 1> StatePVAJ;

class ObvpSolver {
private:
    double rho_;
    double vel_max_ = DBL_MAX, acc_max_ = DBL_MAX;
    double t_star_;
    double cost_star_;
    int test_cnt_ = 0;
public:
    ObvpSolver(int rho) : rho_(rho) {};

    ObvpSolver(int rho, double vel_max, double acc_max) : rho_(rho), vel_max_(vel_max), acc_max_(acc_max) {};

    ~ObvpSolver() {};

    inline double GetLastCost() {
        return cost_star_;
    }

    inline double GetLastTstar() {
        return t_star_;
    }

    inline double CalcOptimalDuration(StatePVAJ &start_state, StatePVAJ &end_state) {
        Eigen::Array3d p_0 = start_state.head(3);
        Eigen::Array3d v_0 = start_state.segment(3, 3);
        Eigen::Array3d a_0 = start_state.segment(6, 3);
        Eigen::Array3d j_0 = start_state.tail(3);

        Eigen::Array3d p_f = end_state.head(3);
        Eigen::Array3d v_f = end_state.segment(3, 3);
        Eigen::Array3d a_f = end_state.segment(6, 3);
        Eigen::Array3d j_f = end_state.tail(3);

        Eigen::VectorXd coeffsGradT(9);

        coeffsGradT(0) = rho_;
        coeffsGradT(2) = (-8 * j_0 * j_f - 16 * j_0.square() - 16 * j_f.square()).sum();
        coeffsGradT(3) = (240 * a_f * j_0 - 240 * a_0 * j_f - 480 * a_0 * j_0 + 480 * a_f * j_f).sum();
        coeffsGradT(4) = (5040 * a_0 * a_f - 2880 * j_0 * v_0 - 2160 * j_0 * v_f - 2160 * j_f * v_0 - 2880 * j_f * v_f -
                          3600 * a_0.square() - 3600 * a_f.square()).sum();
        coeffsGradT(5) = (37440 * a_f * v_0 - 37440 * a_0 * v_f - 43200 * a_0 * v_0 + 43200 * a_f * v_f -
                          6720 * j_0 * p_0 + 6720 * j_0 * p_f - 6720 * j_f * p_0 + 6720 * j_f * p_f).sum();
        coeffsGradT(6) = (100800 * a_0 * p_f - 100800 * a_0 * p_0 + 100800 * a_f * p_0 - 100800 * a_f * p_f -
                          244800 * v_0 * v_f - 129600 * v_0.square() - 129600 * v_f.square()).sum();
        coeffsGradT(7) = (604800 * p_f * v_0 - 604800 * p_0 * v_f - 604800 * p_0 * v_0 + 604800 * p_f * v_f).sum();
        coeffsGradT(8) = (1411200 * p_0 * p_f - 705600 * p_0.square() - 705600 * p_f.square()).sum();

        std::set<double> roots = RootFinder::solvePolynomial(coeffsGradT, DBL_EPSILON, DBL_MAX, 1e-6);

        bool result = false;
        double tau = DBL_MAX;
        double cost = DBL_MAX;

        VectorXd coeffsSnapObjective(7);
        coeffsSnapObjective(0) = (8 * j_0 * j_f + 16 * j_0.square() + 16 * j_f.square()).sum();
        coeffsSnapObjective(1) = (240 * a_0 * j_0 + 120 * a_0 * j_f - 120 * a_f * j_0 - 240 * a_f * j_f).sum();
        coeffsSnapObjective(2) = (960 * j_0 * v_0 - 1680 * a_0 * a_f + 720 * j_0 * v_f + 720 * j_f * v_0 +
                                  960 * j_f * v_f + 1200 * a_0.square() + 1200 * a_f.square()).sum();
        coeffsSnapObjective(3) = (10800 * a_0 * v_0 + 9360 * a_0 * v_f - 9360 * a_f * v_0 - 10800 * a_f * v_f +
                                  1680 * j_0 * p_0 - 1680 * j_0 * p_f + 1680 * j_f * p_0 - 1680 * j_f * p_f).sum();
        coeffsSnapObjective(4) = (20160 * a_0 * p_0 - 20160 * a_0 * p_f - 20160 * a_f * p_0 + 20160 * a_f * p_f +
                                  48960 * v_0 * v_f + 25920 * v_0.square() + 25920 * v_f.square()).sum();
        coeffsSnapObjective(5) = (100800 * p_0 * v_0 + 100800 * p_0 * v_f - 100800 * p_f * v_0 -
                                  100800 * p_f * v_f).sum();
        coeffsSnapObjective(6) = (100800 * p_0.square() - 201600 * p_0 * p_f + 100800 * p_f.square()).sum();
        for (const double &root : roots) {
            double t7 = pow(root, 7);
            double current = rho_ * root + RootFinder::polyVal(coeffsSnapObjective, root) / t7;
            if (current < cost) {
                tau = root;
                cost = current;
                result = true;
            }
        }
        cost_star_ = cost;
        t_star_ = tau;
        return tau;
    }

    inline double EvaluateSnapCost(Piece pie_in) {
        DynamicMat boundCond = pie_in.getBoundCond();
        Eigen::Array3d p_0 = boundCond.block<3, 1>(0, 0);
        Eigen::Array3d v_0 = boundCond.block<3, 1>(0, 1);
        Eigen::Array3d a_0 = boundCond.block<3, 1>(0, 2);
        Eigen::Array3d j_0 = boundCond.block<3, 1>(0, 3);

        Eigen::Array3d p_f = boundCond.block<3, 1>(0, 4);
        Eigen::Array3d v_f = boundCond.block<3, 1>(0, 5);
        Eigen::Array3d a_f = boundCond.block<3, 1>(0, 6);
        Eigen::Array3d j_f = boundCond.block<3, 1>(0, 7);
        VectorXd coeffsSnapObjective(7);
        coeffsSnapObjective(0) = (8 * j_0 * j_f + 16 * j_0.square() + 16 * j_f.square()).sum();
        coeffsSnapObjective(1) = (240 * a_0 * j_0 + 120 * a_0 * j_f - 120 * a_f * j_0 - 240 * a_f * j_f).sum();
        coeffsSnapObjective(2) = (960 * j_0 * v_0 - 1680 * a_0 * a_f + 720 * j_0 * v_f + 720 * j_f * v_0 +
                                  960 * j_f * v_f + 1200 * a_0.square() + 1200 * a_f.square()).sum();
        coeffsSnapObjective(3) = (10800 * a_0 * v_0 + 9360 * a_0 * v_f - 9360 * a_f * v_0 - 10800 * a_f * v_f +
                                  1680 * j_0 * p_0 - 1680 * j_0 * p_f + 1680 * j_f * p_0 - 1680 * j_f * p_f).sum();
        coeffsSnapObjective(4) = (20160 * a_0 * p_0 - 20160 * a_0 * p_f - 20160 * a_f * p_0 + 20160 * a_f * p_f +
                                  48960 * v_0 * v_f + 25920 * v_0.square() + 25920 * v_f.square()).sum();
        coeffsSnapObjective(5) = (100800 * p_0 * v_0 + 100800 * p_0 * v_f - 100800 * p_f * v_0 -
                                  100800 * p_f * v_f).sum();
        coeffsSnapObjective(6) = (100800 * p_0.square() - 201600 * p_0 * p_f + 100800 * p_f.square()).sum();
        double t = pie_in.getDuration();
        double t7 = pow(t, 7);
        double current = rho_ * t + RootFinder::polyVal(coeffsSnapObjective, t) / t7;
        return current;
    }

    inline Piece GenFixStateMinSnapOptT(StatePVAJ &start_state, StatePVAJ &end_state) {
        double t = CalcOptimalDuration(start_state, end_state);
        return GenFixStateMinSnap(start_state, end_state, t);
    }

    inline bool CheckPiece(Piece &pie_in) {
        if (pie_in.checkMaxAccRate(acc_max_) && pie_in.checkMaxVelRate(vel_max_)) {
            return true;
        }
        return false;
    }




    inline Piece GenFixStateMinSnapOptTC(StatePVAJ &start_state, StatePVAJ &end_state) {
//#define TIME_CHEAK
#ifdef TIME_CHEAK
        std::chrono::high_resolution_clock::time_point tc_1, tc_2, tc_3, tc_4;
        tc_1 = std::chrono::high_resolution_clock::now();
#endif

        // 1) Make sure the fix state is feasible
        EnforceFeasibility(start_state, end_state);

        // 2) Gen unconstrained primitive as global minimum
        Piece global_min = GenFixStateMinSnapOptT(start_state, end_state);

        // 3) If feasible, return
        if (CheckPiece(global_min)) {
//            printf(" \033[32m Reach global minuimum\033[0m\n");
            return global_min;
        }
        test_cnt_ = 0;

#ifdef TIME_CHEAK
        tc_2 = std::chrono::high_resolution_clock::now();
#endif
        // 4) Unfeasible, try bisection, first find feasible solution
        double test_t = t_star_*2, t_fea= t_star_*2;
        double scale_ratio = 2;
        while (1) {
//            printf("test_t %lf\n", test_t);
            test_t *= scale_ratio;
            Piece check_pie = GenFixStateMinSnap(start_state, end_state, test_t);
            if (CheckPiece(check_pie)) {
                t_fea = test_t;
                break;
            } else {
//                printf("max_v =  %lf, max_a = %lf\n",check_pie.getMaxVelRate(),check_pie.getMaxAccRate());
            }
        }

        int max_it = 5;
        double t_l = t_star_, t_r = t_fea;
        Piece check_pie;
        bool feasible = false;
#ifdef TIME_CHEAK
        tc_3 = std::chrono::high_resolution_clock::now();
#endif
        while (max_it--) {
            double t_now = (t_l + t_r) / 2;
            Piece check_pie = GenFixStateMinSnap(start_state, end_state, t_now);
            if (CheckPiece(check_pie)) {
                feasible =true;
                t_r = t_now;
            } else {
                t_l = t_now;
                feasible = false;
            }
        }


        check_pie = GenFixStateMinSnap(start_state, end_state, t_r);
        t_star_ = t_r;
        cost_star_ = EvaluateSnapCost(check_pie);

#ifdef TIME_CHEAK
        tc_4 = std::chrono::high_resolution_clock::now();
        double dt = std::chrono::duration_cast<std::chrono::duration<double>>(tc_2 - tc_1).count();
        double t_us = (double) dt * 1e6;
        printf("Preprocess time consuming \033[32m %lf us\033[0m\n", t_us);
        dt = std::chrono::duration_cast<std::chrono::duration<double>>(tc_3 - tc_2).count();
        t_us = (double) dt * 1e6;
        printf("Opt time consuming \033[32m %lf us\033[0m\n", t_us);
        dt = std::chrono::duration_cast<std::chrono::duration<double>>(tc_4 - tc_3).count();
        t_us = (double) dt * 1e6;
        printf("Bisection time consuming \033[32m %lf us\033[0m\n", t_us);
#endif
        return check_pie;
    }

    inline void EnforceFeasibility(StatePVAJ &start_state, StatePVAJ &end_state) {
        if (start_state.segment(3, 3).norm() > vel_max_) {
            start_state.segment(3, 3) = start_state.segment(3, 3).normalized() * vel_max_;
        }

        if (start_state.segment(6, 3).norm() > acc_max_) {
            start_state.segment(6, 3) = start_state.segment(3, 3).normalized() * acc_max_;
        }

        if (end_state.segment(3, 3).norm() > vel_max_) {
            end_state.segment(3, 3) = end_state.segment(3, 3).normalized() * vel_max_;
        }

        if (end_state.segment(6, 3).norm() > acc_max_) {
            end_state.segment(6, 3) = end_state.segment(3, 3).normalized() * acc_max_;
        }
    }

    inline Piece GenFixStateMinSnap(StatePVAJ &start_state, StatePVAJ &end_state, double t) {
        test_cnt_++;
        Eigen::Array3d p_0 = start_state.head(3);
        Eigen::Array3d v_0 = start_state.segment(3, 3);
        Eigen::Array3d a_0 = start_state.segment(6, 3);
        Eigen::Array3d j_0 = start_state.tail(3);

        Eigen::Array3d p_f = end_state.head(3);
        Eigen::Array3d v_f = end_state.segment(3, 3);
        Eigen::Array3d a_f = end_state.segment(6, 3);
        Eigen::Array3d j_f = end_state.tail(3);
        double tvec[8];
        DynamicMat out_mat(3, 8);
        tvec[0] = 1;
        for (size_t i = 1; i < 8; i++) {
            tvec[i] = pow(t, i);
        }
        double J = 0.0;
        for (size_t i = 0; i < 3; i++) {
            Eigen::Vector4d delta_pvaj = Eigen::Vector4d(
                    p_f[i] - p_0[i] - v_0[i] * t - 1.0 / 2.0f * a_0[i] * tvec[2] -
                    1.0 / 6.0f * j_0[i] * tvec[3],
                    v_f[i] - v_0[i] - a_0[i] * t - 1.0 / 2.0f * j_0[i] * tvec[2],
                    a_f[i] - a_0[i] - j_0[i] * t,
                    j_f[i] - j_0[i]);

            Eigen::MatrixXd tmat(4, 4);
            tmat << -33600, 16800 * tvec[1], -3360 * tvec[2], 280 * tvec[3],
                    25200 * tvec[1], -12240 * tvec[2], 2340 * tvec[3], -180 * tvec[4],
                    -10080 * tvec[2], 4680 * tvec[3], -840 * tvec[4], 60 * tvec[5],
                    840 * tvec[3], -360 * tvec[4], 60 * tvec[5], -4 * tvec[6];

            Eigen::Vector4d abyw = tmat * delta_pvaj / tvec[7];
            double aa = abyw[0], bb = abyw[1], yy = abyw[2], ww = abyw[3];
            out_mat.row(i)[0] = aa / 1680.0f;
            out_mat.row(i)[1] = bb / 360.0f;
            out_mat.row(i)[2] = yy / 120.0f;
            out_mat.row(i)[3] = ww / 24.0f;
            out_mat.row(i)[4] = j_0[i] / 6.0f;
            out_mat.row(i)[5] = a_0[i] / 2.0f;
            out_mat.row(i)[6] = v_0[i];
            out_mat.row(i)[7] = p_0[i];
            J += (aa * aa * tvec[7]) / 28.0 + (aa * bb * tvec[6]) / 12.0 + aa * yy * tvec[5] / 5.0 +
                 aa * ww * tvec[4] / 4.0
                 + (bb * bb * tvec[5]) / 20.0 + bb * yy * tvec[4] / 4.0 + bb * ww * tvec[3] / 3.0 +
                 yy * yy * tvec[3] / 3.0 +
                 yy * ww * tvec[2] + ww * ww * t;
//            J += (aa*aa*tvec[9])/576.0 + (aa*bb*tvec[8])/192.0 + aa*yy*tvec[7]/56.0 +
//                    aa*j_0[i]*tvec[5]/20.0 + aa*ww*tvec[6]/24.0 + (bb*bb*tvec[7])/252.0 +
//                    bb*yy*tvec[6]/36.0+ bb*j_0[i]*tvec[4]/12.0+  bb*ww*tvec[5]/15.0 +
//                    yy*yy*tvec[5]/20.0 +yy*j_0[i]*tvec[3]/3.0+ +yy*ww*tvec[4]/4.0+
//                    j_0[i]*j_0[i]*t +  j_0[i]*ww*tvec[2]  + ww*ww*tvec[3]/3.0;
        }

        Piece out_pie(t, out_mat);
        return out_pie;
    }


public:

};


#endif //SRC_OBVP_SOLVER_HPP
