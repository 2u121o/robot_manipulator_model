#ifndef INVERSE_KINEMATIC_H
#define INVERSE_KINEMATIC_H

#include "kinematic_model.h"

class InverseKinematic{

    public:

        InverseKinematic();

        void solveIk(Eigen::VectorXd &solution);

        ~InverseKinematic();

        void setDesiredPos(Eigen::Vector3d &desired_pos);
        void setQk(Eigen::VectorXd &q_k);

    private:

        const int DOFS = 6;
        const double EPSILON_ERROR = 0.001;

        Eigen::Vector3d desired_pos_;

        Eigen::VectorXd q_k_plus_one_;
        Eigen::VectorXd q_k_;

        KinematicModel kinematic_model_;

        void initVariables();

};

#endif