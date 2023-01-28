#include "kinematic_model.h"

KinematicModel::KinematicModel(){

    initVariables();
    
}

void KinematicModel::initVariables(){

    q_.resize(DOFS);
    jacobian_.resize(DOFS, DOFS);
    alpha_.resize(DOFS);
    a_.resize(DOFS);
    d_.resize(DOFS);

    q_.setZero();
    cartesina_pos_.setZero();
    jacobian_.setZero();

    alpha_ << M_PI/2, M_PI, M_PI/2, -M_PI/2, M_PI/2, 0.0;
    a_ << 0.0, 0.38, 0.0, 0.0, 0.0, 0.0;
    d_ << 0.22, 0.0, 0.0, 0.42, 0.0, 0.157;
}

void KinematicModel::computeForwardKinematic(int idx_first_link, int idx_last_link){

    R_.setIdentity();
    trans_ << 0,0,0;

    Eigen::Matrix3d R; 
    Eigen::Vector3d trans;

    double theta;

    for(int i=idx_first_link; i<idx_last_link; i++){

         if(i==1 || i==2)
            theta = q_[i]+M_PI/2;
        else
            theta = q_[i];

        R << cos(theta), -sin(theta)*cos(alpha_[i]), sin(theta)*sin(alpha_[i]),
            sin(theta), cos(theta)*cos(alpha_[i]) , -cos(theta)*sin(alpha_[i]),
            0           , sin(alpha_[i])              , cos(alpha_[i]);

        trans << a_[i]*cos(theta), a_[i]*sin(theta), d_[i];

        trans_ = R_*trans + trans_;
        R_ = R_*R;
    }

}

void KinematicModel::setQ(const Eigen::VectorXd &q){
    q_ = q;
}

Eigen::VectorXd KinematicModel::getQ(){
    return q_;
}

Eigen::Vector3d KinematicModel::getTrans(){
    return trans_;
}

Eigen::Matrix3d KinematicModel::getR(){
    return R_;
}


KinematicModel::~KinematicModel(){

}