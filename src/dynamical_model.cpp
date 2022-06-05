#include "dynamical_model.h"

DynamicalModel::DynamicalModel(){

    //masses
    m_.resize(DOFS);
    m_ << 1, 1, 1;

    //inertias
    inertia_.resize(DOFS);
    inertia_.at(0) << 1, 0, 0,
                    0, 1, 0,
                    0, 0, 1;
    inertia_.at(1) << 1, 0, 0,
                    0, 1, 0,
                    0, 0, 1;
    inertia_.at(2) << 1, 0, 0,
                    0, 1, 0,
                    0, 0, 1;

    //cogs
    cog_.resize(DOFS);
    cog_.at(0) << 0.2, -0.01, -0.05;
    cog_.at(1) << -0.1, 0.0, 0.0;
    cog_.at(2) << 0.0, 0.0, 0.001;

}

void DynamicalModel::forwardRecursion(Eigen::VectorXd q, Eigen::VectorXd dq, Eigen::VectorXd ddq){


}
