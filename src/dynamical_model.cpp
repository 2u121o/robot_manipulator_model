#include "dynamical_model.h"

DynamicalModel::DynamicalModel(){

    //##########Dynamic Parameters #########//
    //masses
    m_.resize(DOFS);
    m_ << 7.9637582, 2.6189549, 7.5440000, 2.3488000, 2.3432207, 0.82857161;

    //inertias
    inertia_.resize(DOFS);
    inertia_.at(0) << 0.047361000, 0, 0,
                      0, 0.027948983, 0,
                      0, 0, 0.040138000;
    inertia_.at(1) << 0.89171709, 0, 0,
                      0, 1.7625421, 0,
                      0, 0, 2.5169801;
    inertia_.at(2) << 0.66028304, 0, 0,
                      0, 0.28351611, 0,
                      0, 0, 0.93963656;

    inertia_.at(3) << 0.0, 0, 0,
                      0, 0.36765057, 0,
                      0, 0, 0.36765057;

    inertia_.at(4) << 0.0, 0, 0,
                      0, 0.0013644056, 0,
                      0, 0, 0.0013644056;

    inertia_.at(5) << 0.014304167, -0.045700676, -0.081859114,
                       -0.045700676, 0.010091649, -0.095639984,
                      -0.081859114, -0.095639984, 0.024395716;



    //cogs expressed wrt the next frame so ready to use in the formulas without transformation
    cog_.resize(DOFS);
    cog_.at(0) << -0.0023983783, -0.051099000, -0.0064075029;
    cog_.at(1) << -0.27199000, -0.050000000, -0.080318000;
    cog_.at(2) << 0.0062422845, 0.030640000, 0.063114058;
    cog_.at(3) << 0.035455313, 0.21296500, -0.10879600;
    cog_.at(4) << 0.0038853698, 0.025798751, 0.08656600;
    cog_.at(5) << 0.0035370346, 0.0021101928, -0.25491905;


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



Eigen::VectorXd DynamicalModel::rnea(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq, const Eigen::Vector3d gravity){

    initializeMatrices(gravity);

    kinematic_model_.setQ(q);

    forwardRecursion(q, dq, ddq);
    backwardRecursion(q);

    return  u_;

}

void DynamicalModel::forwardRecursion(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq){

    Eigen::Vector3d z;
    z << 0, 0, 1;

    Eigen::Matrix3d R;
    Eigen::Matrix3d R_transpose;
    Eigen::Vector3d t;
  
    for(short int i=0; i<DOFS; i++){

        kinematic_model_.setQ(q);
        kinematic_model_.computeForwardKinematic(i,i+1);
        R = kinematic_model_.getR();
        t = kinematic_model_.getTrans();

        R_transpose = R.transpose();

        omega_.at(i+1) = R_transpose*(omega_.at(i) + dq[i]*z);
        d_omega_.at(i+1) = R_transpose*(d_omega_.at(i) + ddq[i]*z + dq[i]*(omega_.at(i).cross(z)));
        a_.at(i+1) = R_transpose*a_.at(i) + d_omega_.at(i+1).cross(R_transpose*t) + omega_.at(i+1).cross(omega_.at(i+1).cross(R_transpose*t));    
        ac_.at(i) = a_.at(i+1) + d_omega_.at(i+1).cross(cog_.at(i)) + omega_.at(i+1).cross(omega_.at(i+1).cross(cog_.at(i)));

    }

}


void DynamicalModel::backwardRecursion(const Eigen::VectorXd &q){


    Eigen::Vector3d z;
    z << 0, 0, 1;

    Eigen::Matrix3d R;
    Eigen::Matrix3d R_transpose;
    Eigen::Vector3d t;
    Eigen::Matrix3d Rp1;

    R.setZero();
    t.setZero();

    //the first transfomation is identity because it transfom
    //from the tip to the environment so it remain like it is
    Rp1 << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;

    for(short int i=DOFS-1; i>=0; i--){

        kinematic_model_.setQ(q);
        kinematic_model_.computeForwardKinematic(i,i+1);
        R = kinematic_model_.getR();
        t = kinematic_model_.getTrans();

        R_transpose = R.transpose();

       f_.at(i) = Rp1*f_.at(i+1) + m_[i]*(ac_.at(i));
       
        tau_.at(i) = Rp1*tau_.at(i+1) + (Rp1*f_.at(i+1)).cross(cog_.at(i)) - f_.at(i).cross(R_transpose*t + cog_.at(i)) +
                    inertia_.at(i)*(d_omega_.at(i+1)) + omega_.at(i+1).cross(inertia_.at(i) * (omega_.at(i+1)));

        u_[i] = tau_.at(i).transpose()*R_transpose*z;
        Rp1 = R;
    }

}

void DynamicalModel::initializeMatrices(const Eigen::Vector3d gravity){

    for(short int i=0; i<DOFS+1; i++){
        omega_.at(i).setZero();
        d_omega_.at(i).setZero();
        a_.at(i).setZero();
        f_.at(i).setZero();
        tau_.at(i).setZero();
        if(i<DOFS) ac_.at(i).setZero();
    }
    a_.at(0) = gravity;

    u_.setZero();
    
}


