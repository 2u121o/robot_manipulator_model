#include "kinematic_model.h"

KinematicModel::KinematicModel(){

    initVariables();
    std::cout << "Kinematic model constructed! " << std::endl;
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

    // alpha_ << M_PI/2, M_PI, M_PI/2, -M_PI/2, M_PI/2, 0.0;
    // a_ << 0.0, 0.38, 0.0, 0.0, 0.0, 0.0;
    // d_ << 0.22, 0.0, 0.0, 0.42, 0.0, 0.157;

    //UR5 from https://www.universal-robots.com/articles/ur/application-installation/dh-parameters-for-calculations-of-kinematics-and-dynamics/
    //change also in computeForwardKinematic() if change robot
    alpha_ << M_PI/2, 0, 0, M_PI/2, -M_PI/2, 0.0;
    a_ << 0.0, -0.425, -0.39225, 0.0, 0.0, 0.0;
    d_ << 0.089159, 0.0, 0.0, 0.10915, 0.09465, 0.0823;
}

void KinematicModel::computeForwardKinematic(int idx_first_link, int idx_last_link){
 
    R_.setIdentity();
    trans_.setZero();

    Eigen::Matrix3d R; 
    R.setZero();
    Eigen::Vector3d trans;
    trans.setZero();

    double theta;

    for(int i=idx_first_link; i<idx_last_link; i++){

        //  if(i==1 || i==2)
        //     theta = q_[i]+M_PI/2;
        // else
            theta = q_[i];

        R << cos(theta), -sin(theta)*cos(alpha_[i]), sin(theta)*sin(alpha_[i]),
            sin(theta), cos(theta)*cos(alpha_[i]) , -cos(theta)*sin(alpha_[i]),
            0           , sin(alpha_[i])              , cos(alpha_[i]);

        trans << a_[i]*cos(theta), a_[i]*sin(theta), d_[i];

        trans_ = R_*trans + trans_;
        R_ = R_*R;
    }

}

void KinematicModel::computeJacobian(){

    Eigen::Vector3d pos_ee_absolute;
    pos_ee_absolute.setZero();
    computeForwardKinematic(0, DOFS);
    pos_ee_absolute = trans_;

    Eigen::Vector3d pos_ee_relative;
    pos_ee_relative.setZero();

    Eigen::Vector3d z_i_minus_one;
    z_i_minus_one.setZero();

    Eigen::Vector3d z_i_minus_one_const;
    z_i_minus_one_const << 0,0,1;

    for(int i=0; i<DOFS; i++){
        computeForwardKinematic(0, i); //lerrore e probabilmente qui verifica questa parte 
        z_i_minus_one = R_*z_i_minus_one_const;
        pos_ee_relative = pos_ee_absolute - trans_;
        jacobian_.block(0,i,3,1) = z_i_minus_one.cross(pos_ee_relative);
        jacobian_.block(3,i,3,1) = z_i_minus_one;
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

Eigen::MatrixXd KinematicModel::getJacobian(){
    return jacobian_;
}

KinematicModel::~KinematicModel(){

}