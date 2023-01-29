#include "inverse_kinematic.h"

InverseKinematic::InverseKinematic(){
    initVariables();
}

void InverseKinematic::solveIk(Eigen::VectorXd &solution){

    Eigen::MatrixXd jacobian;
    kinematic_model_.setQ(q_k_);
    kinematic_model_.computeJacobian();
    jacobian = kinematic_model_.getJacobian();

    int first_link = 0;
    int last_link = 6;

    kinematic_model_.computeForwardKinematic(first_link, last_link);
    Eigen::Vector3d pos_error = desired_pos_ - kinematic_model_.getTrans();

    while(pos_error.norm()>EPSILON_ERROR){
       
        Eigen::MatrixXd analitic_jacobian = jacobian.block(0,0,3,DOFS);
        Eigen::MatrixXd pinv_jacobian = analitic_jacobian.completeOrthogonalDecomposition().pseudoInverse();
        q_k_plus_one_ = q_k_ + pinv_jacobian*pos_error;
        kinematic_model_.setQ(q_k_plus_one_);
        kinematic_model_.computeForwardKinematic(first_link, last_link);
        pos_error = desired_pos_ - kinematic_model_.getTrans();

        std::cout << "pos_error.norm --> " << pos_error.norm() << std::endl;
    }

    solution = q_k_;

}

void InverseKinematic::initVariables(){

    q_k_plus_one_.resize(DOFS);
    q_k_.resize(DOFS);

    desired_pos_.setZero();
    q_k_plus_one_.setZero();
    q_k_.setZero();
}

void InverseKinematic::setDesiredPos(Eigen::Vector3d &desired_pos){
    desired_pos_ = desired_pos;
}

void InverseKinematic::setQk(Eigen::VectorXd &q_k){
    q_k_ = q_k;
}

InverseKinematic::~InverseKinematic(){

}