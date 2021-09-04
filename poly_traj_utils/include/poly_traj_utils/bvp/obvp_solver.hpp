//
// Created by yunfan on 2021/3/19.
//

#ifndef SRC_OBVP_SOLVER_HPP
#define SRC_OBVP_SOLVER_HPP

#include "../root_finder.hpp"
#include "../poly_data_structure.hpp"
#include "Eigen/Dense"
#include "vector"
#include "../scope_timer.hpp"
#include "fmt/color.h"

using namespace Eigen;
using namespace std;
using namespace fmt;

#define FULLFIX 1
#define FIXPV 2
#define FIXP 3

class ObvpSolver {

private:
    double rho_;
    double vel_max_ = DBL_MAX, acc_max_ = DBL_MAX;
    double t_star_;
    double cost_star_, last_cost_;
    int test_cnt_ = 0;

    inline void PrintState(StatePVAJ &state_) {
        fmt::print(fg(fmt::color::light_yellow), "-- [State]:\n");
        fmt::print("    position:({0:.2f},{1:.2f},{2:.2f})\n", state_.col(0)[0], state_.col(0)[1], state_.col(0)[2]);
        fmt::print("    velocity:({0:.2f},{1:.2f},{2:.2f})\n", state_.col(1)[0], state_.col(1)[1], state_.col(1)[2]);
        fmt::print("    acceleration:({0:.2f},{1:.2f},{2:.2f})\n", state_.col(2)[0], state_.col(2)[1], state_.col(2)[2]);
        fmt::print("    jerk:({0:.2f},{1:.2f},{2:.2f})\n", state_.col(3)[0], state_.col(3)[1], state_.col(3)[2]);
    }

    inline double CalcOptimalDuration(StatePVAJ &start_state, StatePVAJ &end_state, int type, int order = 7) {

        if(order == 7){
            Eigen::Array3d p_0 = start_state.col(0);
            Eigen::Array3d v_0 = start_state.col(1);
            Eigen::Array3d a_0 = start_state.col(2);
            Eigen::Array3d j_0 = start_state.col(3);

            Eigen::Array3d p_f = end_state.col(0);
            Eigen::Array3d v_f = end_state.col(1);
            Eigen::Array3d a_f = end_state.col(2);
            Eigen::Array3d j_f = end_state.col(3);
            Eigen::VectorXd coeffsGradT(9);

            if (type == FULLFIX) {
                coeffsGradT(0) = rho_;
                coeffsGradT(2) = (-8 * j_0 * j_f - 16 * j_0.square() - 16 * j_f.square()).sum();
                coeffsGradT(3) = (240 * a_f * j_0 - 240 * a_0 * j_f - 480 * a_0 * j_0 + 480 * a_f * j_f).sum();
                coeffsGradT(4) = (5040 * a_0 * a_f - 2880 * j_0 * v_0 - 2160 * j_0 * v_f - 2160 * j_f * v_0 -
                        2880 * j_f * v_f -
                        3600 * a_0.square() - 3600 * a_f.square()).sum();
                coeffsGradT(5) = (37440 * a_f * v_0 - 37440 * a_0 * v_f - 43200 * a_0 * v_0 + 43200 * a_f * v_f -
                        6720 * j_0 * p_0 + 6720 * j_0 * p_f - 6720 * j_f * p_0 + 6720 * j_f * p_f).sum();
                coeffsGradT(6) = (100800 * a_0 * p_f - 100800 * a_0 * p_0 + 100800 * a_f * p_0 - 100800 * a_f * p_f -
                        244800 * v_0 * v_f - 129600 * v_0.square() - 129600 * v_f.square()).sum();
                coeffsGradT(7) = (604800 * p_f * v_0 - 604800 * p_0 * v_f - 604800 * p_0 * v_0 + 604800 * p_f * v_f).sum();
                coeffsGradT(8) = (1411200 * p_0 * p_f - 705600 * p_0.square() - 705600 * p_f.square()).sum();
            } else if (type == FIXPV) {

                coeffsGradT(0) = rho_;
                coeffsGradT(2) = (-12.0 * j_0.square()).sum();
                coeffsGradT(3) = (-264.0 * a_0 * j_0).sum();
                coeffsGradT(4) = (-1404 * a_0.square() - 1152 * j_0 * v_0 - 360 * j_0 * v_f).sum();
                coeffsGradT(5) = (2016 * j_0 * p_f - 4320 * a_0 * v_f - 2016 * j_0 * p_0 - 11808 * a_0 * v_0).sum();
                coeffsGradT(6) = (-23760 * v_0.square() - 18000 * v_0 * v_f - 3600 * v_f.square() -
                        20160 * a_0 * p_0 + 20160 * a_0 * p_f).sum();
                coeffsGradT(7) = (78624 * p_f * v_0 - 30240 * p_0 * v_f - 78624 * p_0 * v_0 +
                        30240 * p_f * v_f).sum();
                coeffsGradT(8) = (63504 * p_0 * p_f - 127008 * p_0.square() - 63504 * p_f.square()).sum();
            } else if (type == FIXP) {
                coeffsGradT(0) = rho_;
                coeffsGradT(2) = (-7.0 * j_0.square()).sum();
                coeffsGradT(3) = (-84.0 * a_0 * j_0).sum();
                coeffsGradT(4) = (-189 * a_0.square() - 252 * j_0 * v_0).sum();
                coeffsGradT(5) = (336 * j_0 * p_f - 336 * j_0 * p_0 - 1008 * a_0 * v_0).sum();
                coeffsGradT(6) = (-1260 * v_0.square() - 1260 * a_0 * p_0 + 1260 * a_0 * p_f).sum();
                coeffsGradT(7) = (3024 * p_f * v_0 - 3024 * p_0 * v_0).sum();
                coeffsGradT(8) = (3528 * p_0 * p_f - 1764 * p_0.square() - 1764 * p_f.square()).sum();
            } else {
                fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "Error, undefined bvp type!\n");
            }


            std::set<double> roots = RootFinder::solvePolynomial(coeffsGradT, DBL_EPSILON, DBL_MAX, 1e-6);
            if (roots.empty()) {
                return -1;
            }
            bool result = false;
            double tau = DBL_MAX;
            double cost = DBL_MAX;

            VectorXd coeffsSnapObjective(7);
            if ((type == FULLFIX)) {
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
            } else if ((type == FIXPV)) {

                coeffsSnapObjective(0) = (12 * j_0.square()).sum();
                coeffsSnapObjective(1) = (132 * a_0 * j_0).sum();
                coeffsSnapObjective(2) = (468 * a_0.square() + 384 * j_0 * v_0 + 120 * j_0 * v_f).sum();
                coeffsSnapObjective(3) = (2952 * a_0 * v_0 + 1080 * a_0 * v_f + 504 * j_0 * p_0 -
                        504 * j_0 * p_f).sum();
                coeffsSnapObjective(4) = (4752 * v_0.square() + 3600 * v_0 * v_f + 720 * v_f.square() +
                        4032 * a_0 * p_0 - 4032 * a_0 * p_f).sum();
                coeffsSnapObjective(5) = (13104 * p_0 * v_0 + 5040 * p_0 * v_f - 13104 * p_f * v_0 -
                        5040 * p_f * v_f).sum();
                coeffsSnapObjective(6) = (9072 * p_0.square() - 18144 * p_0 * p_f + 9072 * p_f.square()).sum();
            } else if (type == FIXP) {
                coeffsSnapObjective(0) = (7 * j_0.square()).sum();
                coeffsSnapObjective(1) = (42 * a_0 * j_0).sum();
                coeffsSnapObjective(2) = (63 * a_0.square() + 84 * j_0 * v_0).sum();
                coeffsSnapObjective(3) = (252 * a_0 * v_0 + 84 * j_0 * p_0 - 84 * j_0 * p_f).sum();
                coeffsSnapObjective(4) = (252 * v_0.square() + 252 * a_0 * p_0 - 252 * a_0 * p_f).sum();
                coeffsSnapObjective(5) = (504 * p_0 * v_0 - 504 * p_f * v_0).sum();
                coeffsSnapObjective(6) = (252 * p_0.square() - 504 * p_0 * p_f + 252 * p_f.square()).sum();
            }

            for (const double &root : roots) {
                //            fmt::print(fg(fmt::color::light_pink), "Root: {}\n", root);
                double t7 = pow(root, 7);
                double current = rho_ * root + RootFinder::polyVal(coeffsSnapObjective, root) / t7;
                if (current < cost) {
                    tau = root;
                    cost = current;
                    result = true;
                }
            }
            cost_star_ = cost;
            last_cost_ = cost;
            t_star_ = tau;
            return tau;
        }
        else if(order == 5){

            if(type = FULLFIX){

                Eigen::Vector3d x0_012 = start_state.col(0), x0_345 = start_state.col(1), x0_678 = start_state.col(2);
                Eigen::Vector3d x1_012 = end_state.col(0), x1_345 = end_state.col(1), x1_678 = end_state.col(2);
                Eigen::Vector3d x0_diff_x1_012 = x0_012 - x1_012, x0_plus_x1_345 = 5*x0_345 + 3*x1_345;
                double t1 = 10.0/rho_*x0_678.dot(x0_678);
                double t2 = 10.0/rho_*(7*x0_345.dot(x0_678) + 3*x0_678.dot(x1_345));
                double t3 = 10.0/rho_*(8*x0_345.dot(x0_345) + 3*x1_345.dot(x1_345) + 5*x0_678.dot(x0_diff_x1_012) + 9*x0_345.dot(x1_345));
                double t4 = 10.0/rho_*x0_diff_x1_012.dot(x0_plus_x1_345);
                double t5 = 10.0/rho_*x0_diff_x1_012.dot(x0_diff_x1_012);

                VectorXd p(7);
                p[0] = 1.0;
                p[1] = 0.0;
                p[2] = - 8*t1;
                p[3] = - 16*t2;
                p[4] = - 48*t3;
                p[5] = - 320*t4;
                p[6] = - 1600*t5;
                std::set<double> roots = RootFinder::solvePolynomial(p, DBL_EPSILON, DBL_MAX, 1e-6);

                bool result = false;
                double tau = DBL_MAX;
                double cost = DBL_MAX;

                for (const double& root : roots)
                {
                    double root2 = root*root;
                    double root3 = root2*root;
                    double root4 = root3*root;
                    double root5 = root4*root;
                    double root6 = root5*root;

                    double current = (root6 + 320*t5 + 80*root*t4 + 16*root2*t3 + 8*root3*t2 + 8*root4*t1) / root5;

                    if (current < cost)
                    {
                        tau = root;
                        cost = current;
                        result = true;
                    }
                }
                cost_star_ = cost;
                last_cost_ = cost;
                t_star_ = tau;
                return tau;
            }
            else if( type == FIXPV){

                Eigen::Vector3d x0_012 = start_state.col(0), x0_345 = start_state.col(1), x0_678 = start_state.col(2);
                Eigen::Vector3d x1_012 = end_state.col(0), x1_345 = end_state.col(1), x1_678 = end_state.col(2);
                Eigen::Vector3d x0_diff_x1_012 = x0_012 - x1_012;

                double t4 = 10.0/rho_*x0_678.dot(x0_678);
                double t3 = 10.0/rho_*x0_345.dot(x0_678);
                double t2 = 10.0/rho_*(x0_345.dot(x0_345) + x0_012.dot(x1_678) - x0_678.dot(x1_012));
                double t1 = 10.0/rho_*(x0_012.dot(x0_345) - x0_345.dot(x1_012));
                double t0 = 10.0/rho_*x0_diff_x1_012.dot(x0_diff_x1_012);

                VectorXd p(7);
                p[0] = 1.0;
                p[1] = 0.0;
                p[2] = - 5*t4;
                p[3] = - 40*t3;
                p[4] = - 60*t2;
                p[5] = - 160*t1;
                p[6] = - 100*t0;
                std::set<double> roots = RootFinder::solvePolynomial(p, DBL_EPSILON, DBL_MAX, 1e-6);
                bool result = false;
                double tau = DBL_MAX;
                double cost = DBL_MAX;

                for (const double& root : roots)
                {
                    double root2 = root*root;
                    double root3 = root2*root;
                    double root4 = root3*root;
                    double root5 = root4*root;
                    double root6 = root5*root;

                    double current = (root6 + 20*t0 + 40*root*t1 + 20*root2*t2 + 20*root3*t3 + 5*root4*t4) / root5;
                    if (current < cost)
                    {
                        tau = root;
                        cost = current;
                        result = true;
                    }
                }
                cost_star_ = cost;
                last_cost_ = cost;
                t_star_ = tau;
                return tau;
            }
            else if(type == FIXP){

                Eigen::Vector3d x0_012 = start_state.col(0), x0_345 = start_state.col(1), x0_678 = start_state.col(2);
                Eigen::Vector3d x1_012 = end_state.col(0), x1_345 = end_state.col(1), x1_678 = end_state.col(2);
                Eigen::Vector3d x0_diff_x1_012 = x0_012 - x1_012;

                double t4 = 10.0/rho_*x0_678.dot(x0_678);
                double t3 = 10.0/rho_*x0_345.dot(x0_678);
                double t2 = 10.0/rho_*(x0_345.dot(x0_345) + x0_012.dot(x1_678) - x0_678.dot(x1_012));
                double t1 = 10.0/rho_*(x0_012.dot(x0_345) - x0_345.dot(x1_012));
                double t0 = 10.0/rho_*x0_diff_x1_012.dot(x0_diff_x1_012);

                VectorXd p(7);
                p[0] = 1.0;
                p[1] = 0.0;
                p[2] = - 5*t4;
                p[3] = - 40*t3;
                p[4] = - 60*t2;
                p[5] = - 160*t1;
                p[6] = - 100*t0;
                std::set<double> roots = RootFinder::solvePolynomial(p, DBL_EPSILON, DBL_MAX, 1e-6);
                bool result = false;
                double tau = DBL_MAX;
                double cost = DBL_MAX;

                for (const double& root : roots)
                {
                    double root2 = root*root;
                    double root3 = root2*root;
                    double root4 = root3*root;
                    double root5 = root4*root;
                    double root6 = root5*root;

                    double current = (root6 + 20*t0 + 40*root*t1 + 20*root2*t2 + 20*root3*t3 + 5*root4*t4) / root5;
                    if (current < cost)
                    {
                        tau = root;
                        cost = current;
                        result = true;
                    }
                }

                cost_star_ = cost;
                last_cost_ = cost;
                t_star_ = tau;
                return tau;
            }

        }


    }

    inline void EnforceFeasibility(StatePVAJ &start_state, StatePVAJ &end_state) {
        if (start_state.col(1).norm() > vel_max_) {
            start_state.col(1)= start_state.col(1).normalized() * vel_max_;
        }

        if (start_state.col(2).norm() > acc_max_) {
            start_state.col(2) = start_state.col(2).normalized() * acc_max_;
        }

        if (end_state.col(1).norm() > vel_max_) {
            end_state.col(1)= end_state.col(1).normalized() * vel_max_;
        }

        if (end_state.col(2).norm() > acc_max_) {
            end_state.col(2) = end_state.col(2).normalized() * acc_max_;
        }
    }

    inline bool CheckPiece(Piece &pie_in) {
        if (pie_in.checkMaxAccRate(acc_max_) && pie_in.checkMaxVelRate(vel_max_)) {
            return true;
        }
        return false;
    }

public:
    ObvpSolver(int rho) : rho_(rho) {};

    ObvpSolver(int rho, double vel_max, double acc_max) : rho_(rho), vel_max_(vel_max), acc_max_(acc_max) {

        print(fg(color::yellow_green), " -- [OBVP] param init process\n");
        print(fg(color::dim_gray), " \t\tparam: rho : {}\n", rho);
        print(fg(color::dim_gray), " \t\tparam: max_v : {}\n", vel_max_);
        print(fg(color::dim_gray), " \t\tparam: max_a : {}\n", acc_max_);
    };

    ~ObvpSolver() {};

    typedef shared_ptr<ObvpSolver> Ptr;

    inline double GetLastCostStar() {
        return cost_star_;
    }

    inline double GetLastCost() {
        return last_cost_;
    }

    inline double GetLastTstar() {
        return t_star_;
    }

    inline double EvaluateSnapCost(Piece pie_in, int type = FULLFIX) {
        double t = pie_in.getDuration();
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

        if ((type == FULLFIX)) {
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
        } else if ((type == FIXPV)) {

            coeffsSnapObjective(0) = (12 * j_0.square()).sum();
            coeffsSnapObjective(1) = (132 * a_0 * j_0).sum();
            coeffsSnapObjective(2) = (468 * a_0.square() + 384 * j_0 * v_0 + 120 * j_0 * v_f).sum();
            coeffsSnapObjective(3) = (2952 * a_0 * v_0 + 1080 * a_0 * v_f + 504 * j_0 * p_0 -
                    504 * j_0 * p_f).sum();
            coeffsSnapObjective(4) = (4752 * v_0.square() + 3600 * v_0 * v_f + 720 * v_f.square() +
                    4032 * a_0 * p_0 - 4032 * a_0 * p_f).sum();
            coeffsSnapObjective(5) = (13104 * p_0 * v_0 + 5040 * p_0 * v_f - 13104 * p_f * v_0 -
                    5040 * p_f * v_f).sum();
            coeffsSnapObjective(6) = (9072 * p_0.square() - 18144 * p_0 * p_f + 9072 * p_f.square()).sum();
        } else if (type == FIXP) {
            coeffsSnapObjective(0) = (7 * j_0.square()).sum();
            coeffsSnapObjective(1) = (42 * a_0 * j_0).sum();
            coeffsSnapObjective(2) = (63 * a_0.square() + 84 * j_0 * v_0).sum();
            coeffsSnapObjective(3) = (252 * a_0 * v_0 + 84 * j_0 * p_0 - 84 * j_0 * p_f).sum();
            coeffsSnapObjective(4) = (252 * v_0.square() + 252 * a_0 * p_0 - 252 * a_0 * p_f).sum();
            coeffsSnapObjective(5) = (504 * p_0 * v_0 - 504 * p_f * v_0).sum();
            coeffsSnapObjective(6) = (252 * p_0.square() - 504 * p_0 * p_f + 252 * p_f.square()).sum();
        }

        double t7 = pow(t, 7);
        double current = rho_ * t + RootFinder::polyVal(coeffsSnapObjective, t) / t7;
        return current;
    }

    inline Piece GenFixStateMinJerkOptT(StatePVAJ &start_state, StatePVAJ &end_state, int type = FULLFIX){
        double t = CalcOptimalDuration(start_state, end_state, type, 5);
        if (t < 0 || t > 1000) {
            Piece empt;
            fmt::print(fg(fmt::color::yellow), "Calc Opt T Failed.\n");
            return empt;
        }
        return GenFixStateMinJerk(start_state, end_state, t, type);
    }

    inline Piece GenFixStateMinJerk(StatePVAJ &start_state, StatePVAJ &goal_state, double t, int type = FULLFIX){
        Map<VectorXd> x0_(start_state.data(), 12,1);
        Map<VectorXd> x1_(goal_state.data(), 12,1);
        if(type == FULLFIX){
            DynamicMat coeff_;
            coeff_.resize(3,6);
            double tau_star_ = t;
            double t2 = tau_star_*tau_star_;
            double t3 = tau_star_*t2;
            double t4 = tau_star_*t3;
            double t5 = tau_star_*t4;
            coeff_(0, 0) = -(12*x0_[0] + 6*tau_star_*x0_[3] + t2*x0_[6] - 12*x1_[0] + 6*tau_star_*x1_[3] - t2*x1_[6])/(2*t5);
            coeff_(1, 0) = -(12*x0_[1] + 6*tau_star_*x0_[4] + t2*x0_[7] - 12*x1_[1] + 6*tau_star_*x1_[4] - t2*x1_[7])/(2*t5);
            coeff_(2, 0) = -(12*x0_[2] + 6*tau_star_*x0_[5] + t2*x0_[8] - 12*x1_[2] + 6*tau_star_*x1_[5] - t2*x1_[8])/(2*t5);
            coeff_(0, 1) = -(-30*x0_[0] - 16*tau_star_*x0_[3] - 3*t2*x0_[6] + 30*x1_[0] - 14*tau_star_*x1_[3] + 2*t2*x1_[6])/(2*t4);
            coeff_(1, 1) = -(-30*x0_[1] - 16*tau_star_*x0_[4] - 3*t2*x0_[7] + 30*x1_[1] - 14*tau_star_*x1_[4] + 2*t2*x1_[7])/(2*t4);
            coeff_(2, 1) = -(-30*x0_[2] - 16*tau_star_*x0_[5] - 3*t2*x0_[8] + 30*x1_[2] - 14*tau_star_*x1_[5] + 2*t2*x1_[8])/(2*t4);
            coeff_(0, 2) = -(20*x0_[0] + 12*tau_star_*x0_[3] + 3*t2*x0_[6] - 20*x1_[0] + 8*tau_star_*x1_[3] - t2*x1_[6])/(2*t3);
            coeff_(1, 2) = -(20*x0_[1] + 12*tau_star_*x0_[4] + 3*t2*x0_[7] - 20*x1_[1] + 8*tau_star_*x1_[4] - t2*x1_[7])/(2*t3);
            coeff_(2, 2) = -(20*x0_[2] + 12*tau_star_*x0_[5] + 3*t2*x0_[8] - 20*x1_[2] + 8*tau_star_*x1_[5] - t2*x1_[8])/(2*t3);
            coeff_(0, 3) = x0_[6]/2;
            coeff_(1, 3) = x0_[7]/2;
            coeff_(2, 3) = x0_[8]/2;
            coeff_(0, 4) = x0_[3];
            coeff_(1, 4) = x0_[4];
            coeff_(2, 4) = x0_[5];
            coeff_(0, 5) = x0_[0];
            coeff_(1, 5) = x0_[1];
            coeff_(2, 5) = x0_[2];
            Piece out_pie( t, coeff_);
            out_pie.setCost(last_cost_);
            return out_pie;
        }
        else if(type == FIXPV){
            DynamicMat coeff_;
            coeff_.resize(3,6);
            double tau_star_ = t;
            double t2 = tau_star_*tau_star_;
            double t3 = tau_star_*t2;
            double t4 = tau_star_*t3;
            double t5 = tau_star_*t4;
            coeff_(0, 0) = -(8*x0_[0] + 5*tau_star_*x0_[3] + t2*x0_[6] - 8*x1_[0] + 3*tau_star_*x1_[3])/(3*t5);
            coeff_(1, 0) = -(8*x0_[1] + 5*tau_star_*x0_[4] + t2*x0_[7] - 8*x1_[1] + 3*tau_star_*x1_[4])/(3*t5);
            coeff_(2, 0) = -(8*x0_[2] + 5*tau_star_*x0_[5] + t2*x0_[8] - 8*x1_[2] + 3*tau_star_*x1_[5])/(3*t5);
            coeff_(0, 1) = -(-50*x0_[0] - 32*tau_star_*x0_[3] - 7*t2*x0_[6] + 50*x1_[0] - 18*tau_star_*x1_[3])/(6*t4);
            coeff_(1, 1) = -(-50*x0_[1] - 32*tau_star_*x0_[4] - 7*t2*x0_[7] + 50*x1_[1] - 18*tau_star_*x1_[4])/(6*t4);
            coeff_(2, 1) = -(-50*x0_[2] - 32*tau_star_*x0_[5] - 7*t2*x0_[8] + 50*x1_[2] - 18*tau_star_*x1_[5])/(6*t4);
            coeff_(0, 2) = -(20*x0_[0] + 14*tau_star_*x0_[3] + 4*t2*x0_[6] - 20*x1_[0] + 6*tau_star_*x1_[3])/(3*t3);
            coeff_(1, 2) = -(20*x0_[1] + 14*tau_star_*x0_[4] + 4*t2*x0_[7] - 20*x1_[1] + 6*tau_star_*x1_[4])/(3*t3);
            coeff_(2, 2) = -(20*x0_[2] + 14*tau_star_*x0_[5] + 4*t2*x0_[8] - 20*x1_[2] + 6*tau_star_*x1_[5])/(3*t3);
            coeff_(0, 3) = x0_[6]/2;
            coeff_(1, 3) = x0_[7]/2;
            coeff_(2, 3) = x0_[8]/2;
            coeff_(0, 4) = x0_[3];
            coeff_(1, 4) = x0_[4];
            coeff_(2, 4) = x0_[5];
            coeff_(0, 5) = x0_[0];
            coeff_(1, 5) = x0_[1];
            coeff_(2, 5) = x0_[2];
            Piece out_pie( t, coeff_);
            out_pie.setCost(last_cost_);
            return out_pie;
        }else if(type == FIXP){
            DynamicMat coeff_;
            coeff_.resize(3,6);
            double tau_star_ = t;
            double t2 = tau_star_*tau_star_;
            double t3 = tau_star_*t2;
            double t4 = tau_star_*t3;
            double t5 = tau_star_*t4;
            coeff_(0, 0) = -(2*x0_[0] + 2*tau_star_*x0_[3] + t2*x0_[6] - 2*x1_[0])/(12*t5);
            coeff_(1, 0) = -(2*x0_[1] + 2*tau_star_*x0_[4] + t2*x0_[7] - 2*x1_[1])/(12*t5);
            coeff_(2, 0) = -(2*x0_[2] + 2*tau_star_*x0_[5] + t2*x0_[8] - 2*x1_[2])/(12*t5);
            coeff_(0, 1) = (10*x0_[0] + 10*tau_star_*x0_[3] + 5*t2*x0_[6] - 10*x1_[0])/(12*t4);
            coeff_(1, 1) = (10*x0_[1] + 10*tau_star_*x0_[4] + 5*t2*x0_[7] - 10*x1_[1])/(12*t4);
            coeff_(2, 1) = (10*x0_[2] + 10*tau_star_*x0_[5] + 5*t2*x0_[8] - 10*x1_[2])/(12*t4);
            coeff_(0, 2) = -(10*x0_[0] + 10*tau_star_*x0_[3] + 5*t2*x0_[6] - 10*x1_[0])/(6*t3);
            coeff_(1, 2) = -(10*x0_[1] + 10*tau_star_*x0_[4] + 5*t2*x0_[7] - 10*x1_[1])/(6*t3);
            coeff_(2, 2) = -(10*x0_[2] + 10*tau_star_*x0_[5] + 5*t2*x0_[8] - 10*x1_[2])/(6*t3);
            coeff_(0, 3) = x0_[6]/2;
            coeff_(1, 3) = x0_[7]/2;
            coeff_(2, 3) = x0_[8]/2;
            coeff_(0, 4) = x0_[3];
            coeff_(1, 4) = x0_[4];
            coeff_(2, 4) = x0_[5];
            coeff_(0, 5) = x0_[0];
            coeff_(1, 5) = x0_[1];
            coeff_(2, 5) = x0_[2];

            Piece out_pie( t, coeff_);
            out_pie.setCost(last_cost_);
            return out_pie;
        }

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
        if (global_min.empty()) {
            return global_min;
        }
        // 3) If feasible, return
        if (CheckPiece(global_min)) {
            printf(" \033[32m Reach global minuimum\033[0m\n");
            return global_min;
        }
        test_cnt_ = 0;

#ifdef TIME_CHEAK
        tc_2 = std::chrono::high_resolution_clock::now();
#endif
        // 4) Unfeasible, try bisection, first find feasible solution
        double test_t = t_star_ * 2, t_fea = t_star_ * 2;
        double scale_ratio = 2;
        int maxit = 100;
        while (maxit--) {
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
        if (maxit == 0) {
            Piece emp;
            return emp;
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
                feasible = true;
                t_r = t_now;
            } else {
                t_l = t_now;
                feasible = false;
            }
        }

        if(t_r > 1000){
            Piece pie2;
            return pie2;
        }
        check_pie = GenFixStateMinSnap(start_state, end_state, t_r, true);
        t_star_ = t_r;
        cost_star_ = last_cost_;

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
//        printf("max_v =  %lf, max_a = %lf\n",check_pie.getMaxVelRate(),check_pie.getMaxAccRate());
        return check_pie;
    }

    inline Piece GenFixStateMinSnapOptT(StatePVAJ &start_state, StatePVAJ &end_state, int type = FULLFIX) {
        double t = CalcOptimalDuration(start_state, end_state, type);
        if (t < 0 || t > 1000) {
            Piece empt;
            fmt::print(fg(fmt::color::yellow), "Calc Opt T Failed.\n");
            return empt;
        }
        return GenFixStateMinSnap(start_state, end_state, t, type, false);
    }

    inline Piece GenFixStateMinSnap(StatePVAJ &start_state, StatePVAJ &end_state, double t, int type = FULLFIX,
                                    bool calc_cost = false) {
        calc_cost = true;
        Eigen::Array3d p_0 = start_state.col(0);
        Eigen::Array3d v_0 = start_state.col(1);
        Eigen::Array3d a_0 = start_state.col(2);
        Eigen::Array3d j_0 = start_state.col(3);

        Eigen::Array3d p_f = end_state.col(0);
        Eigen::Array3d v_f = end_state.col(1);
        Eigen::Array3d a_f = end_state.col(2);
        Eigen::Array3d j_f = end_state.col(3);
        double tvec[8];
        DynamicMat out_mat(3, 8);
        tvec[0] = 1;
        for (size_t i = 1; i < 8; i++) {
            tvec[i] = pow(t, i);
        }
        double J = 0.0, aa, bb, yy, ww;
        for (size_t i = 0; i < 3; i++) {
            switch (type) {
                case FULLFIX: {
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
                    aa = abyw[0];
                    bb = abyw[1];
                    yy = abyw[2];
                    ww = abyw[3];
                    break;
                }
                case FIXP: {
                    aa = (42 * a_0[i]) / tvec[5] + (14 * j_0[i]) / tvec[4] + (84 * p_0[i]) / tvec[7] -
                         (84 * p_f[i]) / tvec[7] + (84 * v_0[i]) / tvec[6];
                    bb = (126 * p_f[i]) / tvec[6] - (21 * j_0[i]) / tvec[3] - (126 * p_0[i]) / tvec[6] -
                         (63 * a_0[i]) / tvec[4] - (126 * v_0[i]) / tvec[5];
                    yy = (63 * a_0[i]) / tvec[3] + (21 * j_0[i]) / tvec[2] + (126 * p_0[i]) / tvec[5] -
                         (126 * p_f[i]) / tvec[5] + (126 * v_0[i]) / tvec[4];
                    ww = (42 * p_f[i]) / tvec[4] - (7 * j_0[i]) / t - (42 * p_0[i]) / tvec[4] -
                         (21 * a_0[i]) / tvec[2] - (42 * v_0[i]) / tvec[3];

                    break;
                }
                case FIXPV: {
                     aa = (672 * a_0[i]) / tvec[5] + (84 * j_0[i]) / tvec[4] + (3024 * p_0[i]) / tvec[7] -
                                (3024 * p_f[i]) / tvec[7] + (2184 * v_0[i]) / tvec[6] + (840 * v_f[i]) / tvec[6];
                     bb = (3276 * p_f[i]) / tvec[6] - (96 * j_0[i]) / tvec[3] - (3276 * p_0[i]) / tvec[6] -
                                (738 * a_0[i]) / tvec[4] - (2376 * v_0[i]) / tvec[5] - (900 * v_f[i]) / tvec[5];
                     yy = (468 * a_0[i]) / tvec[3] + (66 * j_0[i]) / tvec[2] + (2016 * p_0[i]) / tvec[5] -
                                (2016 * p_f[i]) / tvec[5] + (1476 * v_0[i]) / tvec[4] + (540 * v_f[i]) / tvec[4];
                     ww = (252 * p_f[i]) / tvec[4] - (12 * j_0[i]) / t - (252 * p_0[i]) / tvec[4] -
                                (66 * a_0[i]) / tvec[2] - (192 * v_0[i]) / tvec[3] - (60 * v_f[i]) / tvec[3];

                    break;
                }
                default: {
                    fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "Error, undefined bvp type!\n");
                }

            }
            out_mat.row(i)[0] = aa / 1680.0f;
            out_mat.row(i)[1] = bb / 360.0f;
            out_mat.row(i)[2] = yy / 120.0f;
            out_mat.row(i)[3] = ww / 24.0f;
            out_mat.row(i)[4] = j_0[i] / 6.0f;
            out_mat.row(i)[5] = a_0[i] / 2.0f;
            out_mat.row(i)[6] = v_0[i];
            out_mat.row(i)[7] = p_0[i];
            if (calc_cost) {
                J += (aa * aa * tvec[7]) / 28.0 + (aa * bb * tvec[6]) / 6.0 + aa * yy * tvec[5] / 5.0 +
                     aa * ww * tvec[4] / 4.0 + (bb * bb * tvec[5]) / 5.0 + bb * yy * tvec[4] / 2.0 +
                     2.0 * bb * ww * tvec[3] / 3.0 +
                     yy * yy * tvec[3] / 3.0 + yy * ww * tvec[2] + ww * ww * t;
            }
        }
        if (calc_cost) {
            last_cost_ = J + rho_ * t;
//            fmt::print(fg(fmt::color::yellow_green) | fmt::emphasis::bold, "Acc j = {}\n", last_cost_);
        }
        Piece out_pie(t, out_mat);
        out_pie.setCost(last_cost_);
        return out_pie;
    }

};


#endif //SRC_OBVP_SOLVER_HPP
