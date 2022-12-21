#include "dynamical_model.h"

DynamicalModel::DynamicalModel(){

    //##########Dynamic Parameters #########//
    //masses
    m_.resize(DOFS);
     m_ << 7.9637582, 2.6189549, 7.5440000, 2.3488000, 2.3432207, 8.2857161;
     //m_ << 7.9637582, 2.6189549, 7.5440000;

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

    inertia_.at(3) << 0.01, 0, 0,
                      0, 0.36765057, 0,
                      0, 0, 0.36765057;

    inertia_.at(4) << 0.01, 0, 0,
                      0, 0.13644056, 0,
                      0, 0, 0.0013644056;

    inertia_.at(5) << 0.014304167, 0, 0,
                      0, 0.010091649, 0,
                      0, 0, 0.024395716;

    // inertia_.at(0) << 0.1, 0, 0,
    //                   0, 0.1, 0,
    //                   0, 0, 0.1;
    // inertia_.at(1) << 0.1, 0, 0,
    //                   0, 0.1, 0,
    //                   0, 0,0.1;
    // inertia_.at(2) << 0.1, 0, 0,
    //                   0, 0.1, 0,
    //                   0, 0, 0.1;

    // inertia_.at(3) << 0.1, 0, 0,
    //                   0, 0.1, 0,
    //                   0, 0, 0.1;

    // inertia_.at(4) << 0.1, 0, 0,
    //                   0, 0.1, 0,
    //                   0, 0, 0.1;

    // inertia_.at(5) << 0.1, 0, 0,
    //                   0, 0.1, 0,
    //                   0, 0, 0.1;


    //cogs expressed wrt the next frame so ready to use in the formulas without transformation
    cog_.resize(DOFS);
    cog_.at(0) << -0.0023983783, -0.051099000, -0.0064075029;
    cog_.at(1) << -0.27199000, -0.050000000, -0.080318000;
    cog_.at(2) << 0.0062422845, 0.030640000, 0.063114058;
    cog_.at(3) << 0.035455313, 0.21296500, -0.10879600;
    cog_.at(4) << 0.0038853698, 0.025798751, 0.08656600;
    cog_.at(5) << 0.0035370346, 0.0021101928, -0.25491905;


    alpha_.resize(DOFS);
    l_.resize(DOFS);
    d_.resize(DOFS);
    alpha_ << M_PI/2, M_PI, M_PI/2, -M_PI/2, M_PI/2, 0.0;
    l_ << 0.0, 0.38, 0.0, 0.0, 0.0, 0.0;
    d_ << 0.22, 0.0, 0.0, 0.42, 0.0, 0.157;
    // alpha_ << M_PI/2, M_PI, M_PI/2;
    // l_ << 0.0, 0.38, 0.0;
    // d_ << 0.22, 0.0, 0.0;

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



Eigen::VectorXd DynamicalModel::rnea(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq){

    initializeMatrices();
    forwardRecursion(q, dq, ddq);
    backwardRecursion(q);

    return  u_;

}

void DynamicalModel::forwardRecursion(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq){

    Eigen::Vector3d z;
    z << 0, 0, 1;

    Eigen::Matrix3d R;
    Eigen::Vector3d t;

    std::vector<double> dh_parameters;
    dh_parameters.resize(4);

   
    for(short int i=0; i<DOFS; i++){

        if(i==1 || i==2)
            dh_parameters.at(0) = q[i]+M_PI/2;
        else
            dh_parameters.at(0) = q[i];
        dh_parameters.at(1) = d_[i];
        dh_parameters.at(2) = l_[i];
        dh_parameters.at(3) = alpha_[i];

       computeRotationTranslation(R, t, dh_parameters);

        omega_.at(i+1) = R.transpose()*(omega_.at(i) + dq[i]*z);
        d_omega_.at(i+1) = R.transpose()*(d_omega_.at(i) + ddq[i]*z + dq[i]*(omega_.at(i).cross(z)));
        // std::cout << d_omega_.at(i+1) << std::endl;
        //  std::cout  << std::endl;

        a_.at(i+1) = R.transpose()*a_.at(i) + d_omega_.at(i+1).cross(R.transpose()*t) + omega_.at(i+1).cross(omega_.at(i+1).cross(R.transpose()*t));
        
        ac_.at(i) = a_.at(i+1) + d_omega_.at(i+1).cross(cog_.at(i)) + omega_.at(i+1).cross(omega_.at(i+1).cross(cog_.at(i)));

       // std::cout << "t --> " << R.transpose()*t << std::endl;
        // std::cout << "omega_.at(i+1) --> " << omega_.at(i+1) << std::endl
        // << "d_omega_.at(i+1) --> " << d_omega_.at(i+1) << std::endl
        // << "a_.at(i+1) --> " << a_.at(i+1) << std::endl
        // << "ac_.at(i) --> " << ac_.at(i) << std::endl;

    }

}


void DynamicalModel::backwardRecursion(const Eigen::VectorXd &q){

    Eigen::Vector3d g;
    g << 0, 0, -9.81;

    Eigen::Vector3d z;
    z << 0, 0, 1;

    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    Eigen::Matrix3d Rp1;

    R.setZero();
    t.setZero();

    std::vector<double> dh_parameters;
    dh_parameters.resize(4);

    //the first transfomation is identity because it transfom
    //from the tip to the environment so it remain like it is
    Rp1 << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;

    std::cout  << tau_.at(0) << std::endl;
    for(short int i=DOFS-1; i>=0; i--){

        if(i==1 || i==2)
            dh_parameters.at(0) = q[i]+M_PI/2;
        else
            dh_parameters.at(0) = q[i];
        dh_parameters.at(1) = d_[i];
        dh_parameters.at(2) = l_[i];
        dh_parameters.at(3) = alpha_[i];

        computeRotationTranslation(R, t, dh_parameters);

       f_.at(i) = Rp1*f_.at(i+1) + m_[i]*(ac_.at(i));
       
       

        tau_.at(i) = Rp1*tau_.at(i+1) + (Rp1*f_.at(i+1)).cross(cog_.at(i)) - f_.at(i).cross(R.transpose()*t + cog_.at(i)) +
                    inertia_.at(i)*(d_omega_.at(i+1)) + omega_.at(i+1).cross(inertia_.at(i) * (omega_.at(i+1)));

        std::cout  << inertia_.at(i) << std::endl;
        std::cout << std::endl;

        u_[i] = tau_.at(i).transpose()*R.transpose()*z;
        Rp1 = R;
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

