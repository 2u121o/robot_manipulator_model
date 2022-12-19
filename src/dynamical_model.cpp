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

    R_.resize(DOFS);
    t_.resize(DOFS);

    initializeMatrices();

}



Eigen::Vector3d DynamicalModel::rnea(const Eigen::Vector3d &q, const Eigen::Vector3d &dq, const Eigen::Vector3d &ddq){

    computeRotationTranslation(q);

    forwardRecursion(q, dq, ddq);
    backwardRecursion(q);

    return  u_;

}

void DynamicalModel::forwardRecursion(const Eigen::Vector3d &q, const Eigen::Vector3d &dq, const Eigen::Vector3d &ddq){

    
    for(short int i=0; i<DOFS; i++){

        double dq_loc = dq[i];
        double ddq_loc = ddq[i];


        Eigen::Matrix3d Ri_transpose = R_[i].transpose();

        omega_.at(i+1) = Ri_transpose*(omega_.at(i) + dq_loc*z_);

        d_omega_.at(i+1) = Ri_transpose*(d_omega_.at(i) + ddq_loc*z_ + dq_loc*(omega_.at(i).cross(z_)));

        a_.at(i+1) = Ri_transpose*a_.at(i) + d_omega_.at(i+1).cross(Ri_transpose*t_[i]) + omega_.at(i+1).cross(omega_.at(i+1).cross(Ri_transpose*t_[i]));

        ac_.at(i) = a_.at(i+1) + d_omega_.at(i+1).cross(cog_.at(i)) + omega_.at(i+1).cross(omega_.at(i+1).cross(cog_.at(i)));

    }

}


void DynamicalModel::backwardRecursion(const Eigen::Vector3d &q){


    Eigen::Matrix3d Rp1;

    //the first transfomation is identity because it transfom
    //from the tip to the environment so it remain like it is
    Rp1 << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;


    for(short int i=DOFS-1; i>=0; i--){

        f_.at(i) = Rp1*f_.at(i+1) + m_[i]*(ac_.at(i));

        tau_.at(i) = Rp1*tau_.at(i+1) + (Rp1*(f_.at(i+1))).cross(cog_.at(i)) - f_.at(i).cross(R_[i].transpose()*t_[i] + cog_.at(i)) +
                    inertia_.at(i)*(d_omega_.at(i+1)) + omega_.at(i+1).cross(inertia_.at(i) * (omega_.at(i+1)));

        u_[i] = tau_.at(i).transpose()*R_[i].transpose()*z_;

    }

}

void DynamicalModel::initializeMatrices(){

    // for(short int i=0; i<DOFS+1; i++){
    //     omega_.at(i).setZero();
    //     d_omega_.at(i).setZero();
    //     a_.at(i).setZero();
    //     f_.at(i).setZero();
    //     tau_.at(i).setZero();
    //     if(i<DOFS) ac_.at(i).setZero();
    // }
    a_.at(0) << 0, 0, 9.81;
    z_ << 0, 0, 1;
    g_ << 0, 0, -9.81;

    //u_.setZero();

}


void DynamicalModel::computeRotationTranslation(const Eigen::Vector3d &q){

    double cq, ca = 0.0;
    double sq, sa = 0.0;
    for(short int i=0; i<DOFS; i++){
        cq = cos(q[i]);
        sq = sin(q[i]);
        ca = cos(alpha_[i]);
        sa = sin(alpha_[i]);
        R_[i] << cq, -sq*cq, sq*sa,
                 sq, cq*ca , -cq*sa,
                 0 , sa    , ca;

        t_[i] << l_[i]*cq, l_[i]*sq, d_[i];
    }
}

