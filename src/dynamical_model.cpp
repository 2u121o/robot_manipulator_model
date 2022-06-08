#include "dynamical_model.h"

DynamicalModel::DynamicalModel(){

    //##########Dynamic Parameters #########//
    //masses
    m_.resize(DOFS);
    m_ << 1, 0.6, 0.2;

    //inertias
    inertia_.resize(DOFS);
    inertia_.at(0) << 0.1, 0, 0,
                      0, 0.1, 0,
                      0, 0, 0.1;
    inertia_.at(1) << 0.1, 0, 0,
                      0, 0.1, 0,
                      0, 0, 0.1;
    inertia_.at(2) << 0.1, 0, 0,
                      0, 0.1, 0,
                      0, 0, 0.1;

    //cogs expressed wrt the next frame so ready to use in the formulas without transformation
    cog_.resize(DOFS);
    cog_.at(0) << 0.01, 0.1, 0;
    cog_.at(1) << -0.2, 0.0, 0.0;
    cog_.at(2) << -0.05, 0.0, 0.01;

    alpha_.resize(DOFS);
    l_.resize(DOFS);
    d_.resize(DOFS);
    alpha_ << M_PI, 0, 0;
    l_ << 0, 0.6, 0.1;
    d_ << 0.3, 0, 0;

    //##########initialization velocities, accelleration#########//
    omega_.resize(DOFS+1);
    d_omega_.resize(DOFS+1);
    a_.resize(DOFS+1);
    ac_.resize(DOFS);


    //##########initialization forces and torques#########//
    f_.resize(DOFS+1);
    tau_.resize(DOFS+1);
    u_.resize(DOFS);

}



Eigen::VectorXd DynamicalModel::rnea(Eigen::VectorXd q, Eigen::VectorXd dq, Eigen::VectorXd ddq){

    initializeMatrices();
    forwardRecursion(q, dq, ddq);
    backwardRecursion(q);

    return  u_;

}

void DynamicalModel::forwardRecursion(Eigen::VectorXd q, Eigen::VectorXd dq, Eigen::VectorXd ddq){

    Eigen::Vector3d z;
    z << 0, 0, 1;

    Eigen::Matrix3d R;
    Eigen::Vector3d t;

    std::vector<double> dh_parameters;
    dh_parameters.resize(4);


    for(short int i=0; i<DOFS; i++){

        dh_parameters.at(0) = q[i];
        dh_parameters.at(1) = d_[i];
        dh_parameters.at(2) = l_[i];
        dh_parameters.at(3) = alpha_[i];


        computeRotationTranslation(R, t, dh_parameters);


        omega_.at(i+1) = R.transpose()*(omega_.at(i) + dq[i]*z);

        d_omega_.at(i+1) = R.transpose()*(d_omega_.at(i) + ddq[i]*z + dq[i]*(omega_.at(i).cross(z)));

        a_.at(i+1) = R.transpose()*a_.at(i) + d_omega_.at(i+1).cross(R.transpose()*t) + omega_.at(i+1).cross(omega_.at(i+1).cross(R.transpose()*t));

        ac_.at(i) = a_.at(i+1) + d_omega_.at(i+1).cross(cog_.at(i)) + omega_.at(i+1).cross(omega_.at(i+1).cross(cog_.at(i)));

    }

}


void DynamicalModel::backwardRecursion(Eigen::VectorXd q){

    Eigen::Vector3d g;
    g << 0, 0, -9.81;

    Eigen::Vector3d z;
    z << 0, 0, 1;

    Eigen::Matrix3d R, Rp1;
    Eigen::Vector3d t;

    std::vector<double> dh_parameters;
    dh_parameters.resize(4);

    //the first transfomation is identity because it transfom
    //from the tip to the environment so it remain like it is
    Rp1 << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;


    for(short int i=DOFS-1; i>=0; i--){


        dh_parameters.at(0) = q[i];
        dh_parameters.at(1) = d_[i];
        dh_parameters.at(2) = l_[i];
        dh_parameters.at(3) = alpha_[i];


        computeRotationTranslation(R, t, dh_parameters);

        f_.at(i) = Rp1*f_.at(i+1) + m_[i]*(ac_.at(i));

        tau_.at(i) = Rp1*tau_.at(i+1) + (Rp1*(f_.at(i+1))).cross(cog_.at(i)) - f_.at(i).cross(R.transpose()*t + cog_.at(i)) +
                    inertia_.at(i)*(d_omega_.at(i+1)) + omega_.at(i+1).cross(inertia_.at(i) * (omega_.at(i+1)));

        u_[i] = tau_.at(i).transpose()*R.transpose()*z;

        computeRotationTranslation(Rp1, t, dh_parameters);


    }

}

void DynamicalModel::initializeMatrices(){

    for(short int i=0; i<DOFS+1; i++){
        omega_.at(i).setZero();
        d_omega_.at(i).setZero();
        a_.at(i).setZero();
        f_.at(i).setZero();
        tau_.at(i).setZero();
        if(i<DOFS) ac_.at(i).setZero();
    }
    a_.at(0) << 0, 0, 9.81;

    u_.setZero();

}


void DynamicalModel::computeRotationTranslation(Eigen::Matrix3d &R, Eigen::Vector3d &t, std::vector<double> dh_params){

    double theta = dh_params[0];
    double d = dh_params[1];
    double a = dh_params[2];
    double alpha = dh_params[3];

    R << cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha),
         sin(theta), cos(theta)*cos(alpha) , -cos(theta)*sin(alpha),
         0           , sin(alpha)              , cos(alpha);

    t << a*cos(theta), a*sin(theta), d;
}

