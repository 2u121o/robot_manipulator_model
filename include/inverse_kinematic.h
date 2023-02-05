#ifndef INVERSE_KINEMATIC_H
#define INVERSE_KINEMATIC_H

#include "kinematic_model.h"

#include <fstream>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

class InverseKinematic{

    public:

        InverseKinematic();

        void solveIk(Eigen::VectorXd &solution);

        ~InverseKinematic();

        void setDesiredPos(Eigen::Vector3d &desired_pos);
        void setQk(Eigen::VectorXd &q_k);

    private:

        const int DOFS = 6;
        const double EPSILON_ERROR = 0.0001;

        Eigen::Vector3d desired_pos_;

        Eigen::VectorXd q_k_plus_one_;
        Eigen::VectorXd q_k_;

        int iteration = 0;

        KinematicModel kinematic_model_;

        std::ofstream fout;

        void initVariables();

        void computeInitialGuess();
        void computeSolution();

};

#endif