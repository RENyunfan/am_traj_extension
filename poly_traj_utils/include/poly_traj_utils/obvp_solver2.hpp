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
using namespace Eigen;
using namespace std;

typedef Eigen::Matrix<double,12,1> StatePVAJ;

class ObvpSolver {
private:
    inline void CalcOptimalDuration(StatePVAJ & start_state, StatePVAJ & end_state){

    }

public:

};


#endif //SRC_OBVP_SOLVER_HPP
