#include "dynamical_model.h"

DynamicalModel::DynamicalModel(){

    //##########Dynamic Parameters #########//
    //masses
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

    alpha_ << M_PI, 0, 0;
    a_ << 0, 0.6, 0.1;
    d_ << 0.3, 0, 0;

    //##########initialization velocities, accelleration#########//
    omega_.resize(DOFS+1);
    d_omega_.resize(DOFS+1);
    a_i_.resize(DOFS+1);
    ac_i_.resize(DOFS);



    //##########initialization forces and torques#########//
    f_.resize(DOFS+1);
    tau_.resize(DOFS+1);



}

void DynamicalModel::initializeMatrices(){

    for(short int i=0; i<DOFS+1; i++){
        omega_.at(i).setZero();
        d_omega_.at(i).setZero();
        a_i_.at(i).setZero();
        f_.at(i).setZero();
        tau_.at(i).setZero();
        if(i<DOFS) ac_i_.at(i).setZero();
    }
    a_i_.at(0) << 0, 0, 9.81;

    u_.setZero();

}

Eigen::Vector3d DynamicalModel::rnea(Eigen::Vector3d q, Eigen::Vector3d dq, Eigen::Vector3d ddq){

    initializeMatrices();
    forwardRecursion(q, dq, ddq);
    backwardRecursion(q);

    return  u_;

}

void DynamicalModel::forwardRecursion(Eigen::Vector3d q, Eigen::Vector3d dq, Eigen::Vector3d ddq){

    Eigen::Vector3d z;
    z << 0, 0, 1;
    Eigen::Matrix3d R;

    double theta_i;
    double alpha_i;

    Eigen::Vector3d t;


    for(short int i=0; i<DOFS; i++){

        theta_i = q[i];
        alpha_i = alpha_[i];
        //is faster to inertia_define directly the transpose
        R << cos(theta_i), -sin(theta_i)*cos(alpha_i), sin(theta_i)*sin(alpha_i),
             sin(theta_i), cos(theta_i)*cos(alpha_i), -cos(theta_i)*sin(alpha_i),
             0           , sin(alpha_i)             , cos(alpha_i);

        t << a_[i]*cos(theta_i), a_[i]*sin(theta_i), d_[i];

        omega_.at(i+1) = R.transpose()*(omega_.at(i) + dq[i]*z);

        d_omega_.at(i+1) = R.transpose()*(d_omega_.at(i) + ddq[i]*z + dq[i]*(omega_.at(i).cross(z)));

        a_i_.at(i+1) = R.transpose()*a_i_.at(i) + d_omega_.at(i+1).cross(R.transpose()*t) + omega_.at(i+1).cross(omega_.at(i+1).cross(R.transpose()*t));

        ac_i_.at(i) = a_i_.at(i+1) + d_omega_.at(i+1).cross(cog_.at(i)) + omega_.at(i+1).cross(omega_.at(i+1).cross(cog_.at(i)));


    }

}


void DynamicalModel::backwardRecursion(Eigen::Vector3d q){
    Eigen::Vector3d g;
    g << 0, 0, -9.81;

    Eigen::Vector3d z;
    z << 0, 0, 1;

    double theta_i;
    double alpha_i;

    Eigen::Matrix3d R, Rp1;
    //the first transfomation is identity because it transfom
    //from the tip to the environment so it remain like it is
    Rp1 << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;

    Eigen::Vector3d t;

    for(short int i=DOFS-1; i>=0; i--){

        theta_i = q[i];
        alpha_i = alpha_[i];


       R << cos(theta_i), -sin(theta_i)*cos(alpha_i), sin(theta_i)*sin(alpha_i),
            sin(theta_i), cos(theta_i)*cos(alpha_i) , -cos(theta_i)*sin(alpha_i),
            0           , sin(alpha_i)              , cos(alpha_i);

       t << a_[i]*cos(theta_i), a_[i]*sin(theta_i), d_[i];

        f_.at(i) = Rp1*f_.at(i+1) + m_[i]*(ac_i_.at(i));

        tau_.at(i) = Rp1*tau_.at(i+1) + (Rp1*(f_.at(i+1))).cross(cog_.at(i)) - f_.at(i).cross(R.transpose()*t + cog_.at(i)) +
                    inertia_.at(i)*(d_omega_.at(i+1)) + omega_.at(i+1).cross(inertia_.at(i) * (omega_.at(i+1)));

        u_[i] = tau_.at(i).transpose()*R.transpose()*z;

        Rp1 << cos(theta_i), -sin(theta_i)*cos(alpha_i), sin(theta_i)*sin(alpha_i),
                sin(theta_i), cos(theta_i)*cos(alpha_i) , -cos(theta_i)*sin(alpha_i),
                0           , sin(alpha_i)              , cos(alpha_i);


    }

}


