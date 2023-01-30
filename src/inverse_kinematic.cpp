#include "inverse_kinematic.h"

InverseKinematic::InverseKinematic(){
    initVariables();
}

void InverseKinematic::solveIk(Eigen::VectorXd &solution){

    
    kinematic_model_.setQ(q_k_);
    kinematic_model_.computeJacobian();
        
    fout.open("ik_solution.dat");
    if(fout) fout << "#iteration position_error_norm" << std::endl;
    computeInitialGuess();
    std::cout << "q_k after IG --> " << q_k_ << std::endl;
    computeSolution();
    fout.close();
    solution = q_k_;

}

void InverseKinematic::computeInitialGuess(){

    Eigen::MatrixXd jacobian;
    jacobian = kinematic_model_.getJacobian();
    Eigen::VectorXd q(DOFS);
    q = q_k_;

    int first_link = 0;
    int last_link = 6;

    kinematic_model_.computeForwardKinematic(first_link, last_link);
    Eigen::Vector3d pos_error = desired_pos_ - kinematic_model_.getTrans();

    double dpos_error_dt = 0;
    double pos_error_norm_prec = 0;

    bool is_ig_found = false;

    while(pos_error.norm()>EPSILON_ERROR && iteration < 1000 && !is_ig_found){

        Eigen::MatrixXd analitic_jacobian = jacobian.block(0,0,3,DOFS);
        Eigen::VectorXd dq = 0.8*analitic_jacobian.transpose()*pos_error;
        q += 0.1*dq;
        q_k_ = q_k_ + q;
        kinematic_model_.setQ(q);
        kinematic_model_.computeForwardKinematic(first_link, last_link);
        kinematic_model_.computeJacobian();
        jacobian = kinematic_model_.getJacobian();
        pos_error = desired_pos_ - kinematic_model_.getTrans();

        if(fout)
            fout << iteration << " " <<  pos_error.norm() << " ";
       
        iteration++;
        dpos_error_dt = (pos_error.norm()-pos_error_norm_prec)/0.001;
        pos_error_norm_prec = pos_error.norm();

        if(std::abs(dpos_error_dt) < 0.0001){
            is_ig_found = true;
        }
        
        if(fout)
            fout << dpos_error_dt << std::endl;
          
    }

    std::cout << "q_k in IG --> " << q_k_ << std::endl;

}

void InverseKinematic::computeSolution(){
    std::cout << "q_k in compute solution --> " << q_k_ << std::endl;
    for(int i=0; i<DOFS; i++){
        q_k_[i] = std::atan2(sin(q_k_[i]),cos(q_k_[i]));
        if(q_k_[i]<0) q_k_[i] += 2*M_PI;
    }
    Eigen::MatrixXd jacobian;
    kinematic_model_.setQ(q_k_);
    jacobian = kinematic_model_.getJacobian();
    Eigen::VectorXd q(DOFS);
    q = q_k_;

    int first_link = 0;
    int last_link = 6;

    kinematic_model_.computeForwardKinematic(first_link, last_link);
    Eigen::Vector3d pos_error = desired_pos_ - kinematic_model_.getTrans();

    while(pos_error.norm()>EPSILON_ERROR && iteration < 6500 ){

        Eigen::MatrixXd analitic_jacobian = jacobian.block(0,0,3,DOFS);
        Eigen::MatrixXd pinv_jacobian = analitic_jacobian.completeOrthogonalDecomposition().pseudoInverse();
        Eigen::VectorXd dq = pinv_jacobian*pos_error;
        q_k_ +=  0.5*dq;
        kinematic_model_.setQ(q_k_);
        kinematic_model_.computeForwardKinematic(first_link, last_link);
        kinematic_model_.computeJacobian();
        jacobian = kinematic_model_.getJacobian();
        pos_error = desired_pos_ - kinematic_model_.getTrans();

        if(fout)
            fout << iteration << " " <<  pos_error.norm() << std::endl;
       
        iteration++;
        
        // if(fout)
        //     fout << 0 << std::endl;
          
    }
}

// void InverseKinematic::solveIk(Eigen::VectorXd &solution){

//     Eigen::MatrixXd jacobian;
//     kinematic_model_.setQ(q_k_);
//     kinematic_model_.computeJacobian();
//     jacobian = kinematic_model_.getJacobian();

//     int first_link = 0;
//     int last_link = 6;

//     kinematic_model_.computeForwardKinematic(first_link, last_link);
//     Eigen::Vector3d pos_error = desired_pos_ - kinematic_model_.getTrans();

//     std::ofstream fout;
//     fout.open("ik_solution.dat");
//     int iteration = 0;
//     if(fout)
//             fout << "#iteration position_error_norm" << std::endl;
//     while(pos_error.norm()>EPSILON_ERROR && iteration < 200){
       
//         Eigen::MatrixXd analitic_jacobian = jacobian.block(0,0,3,DOFS);
//         //Eigen::MatrixXd pinv_jacobian = analitic_jacobian.completeOrthogonalDecomposition().pseudoInverse();
//         Eigen::MatrixXd pinv_jacobian = analitic_jacobian.transpose();
//         q_k_plus_one_ = q_k_ + 0.2*pinv_jacobian*pos_error;
//         kinematic_model_.setQ(q_k_plus_one_);
//         kinematic_model_.computeForwardKinematic(first_link, last_link);
//         pos_error = desired_pos_ - kinematic_model_.getTrans();
//         q_k_ = q_k_plus_one_;
        
//         if(fout)
//             fout << iteration << " " <<  pos_error.norm() << std::endl;

//         iteration++;
//     }
//     fout.close();
//     solution = q_k_;

// }

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