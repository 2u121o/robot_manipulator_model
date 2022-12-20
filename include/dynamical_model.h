#ifndef DYNAMICAL_MODEL_H
#define DYNAMICAL_MODEL_H

#include <iostream>
#include <vector>
#include <math.h>

#include <Eigen/Dense>

class DynamicalModel{

public:
    DynamicalModel();

    Eigen::VectorXd rnea(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq);

private:
    const int DOFS = 6;
    Eigen::VectorXd m_;

    std::vector<Eigen::Matrix3d> inertia_; //3X3N tensor contains inertia matrix for each link
    std::vector<Eigen::Vector3d> cog_;     //3XN matrix contains cog vector for each link

    //parameters of the DH table
    Eigen::VectorXd alpha_;
    Eigen::VectorXd l_;
    Eigen::VectorXd d_;

    std::vector<Eigen::Vector3d> omega_;    //angular velocity
    std::vector<Eigen::Vector3d> d_omega_;  //angular acceleration
    std::vector<Eigen::Vector3d> a_;      //linear acceleration at frame i in frame i
    std::vector<Eigen::Vector3d> ac_;     //linear acceleration of the cog in frame i

    std::vector<Eigen::Vector3d> f_;        //linear forces applied at the origin of the RF
    std::vector<Eigen::Vector3d> tau_;      //tourques

    Eigen::VectorXd u_;                     //tau_ after projection

    void initializeMatrices();

    void forwardRecursion(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq);

    void backwardRecursion(const Eigen::VectorXd &q);

    //in the dh_parameters the values are in the following order
    //theta d a alpha
    void computeRotationTranslation(Eigen::Matrix3d &R, Eigen::Vector3d &t, std::vector<double> dh_params);



};

#endif
