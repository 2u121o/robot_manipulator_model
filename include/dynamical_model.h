#ifndef DYNAMICAL_MODEL_H
#define DYNAMICAL_MODEL_H

#include <iostream>
#include <vector>
#include <math.h>

#include <Eigen/Dense>

class DynamicalModel{

public:
    DynamicalModel();

    Eigen::Vector3d rnea( Eigen::Vector3d q, Eigen::Vector3d dq, Eigen::Vector3d ddq);

private:
    const int DOFS = 3;
    Eigen::Vector3d m_;

    std::vector<Eigen::Matrix3d> inertia_; //3X3N tensor contains inertia matrix for each link
    std::vector<Eigen::Vector3d> cog_;     //3XN matrix contains cog vector for each link

    //parameters of the DH table
    Eigen::Vector3d alpha_;
    Eigen::Vector3d a_;
    Eigen::Vector3d d_;

    std::vector<Eigen::Vector3d> omega_;
    std::vector<Eigen::Vector3d> d_omega_;
    std::vector<Eigen::Vector3d> a_i_;
    std::vector<Eigen::Vector3d> ac_i_;

    std::vector<Eigen::Vector3d> f_;
    std::vector<Eigen::Vector3d> tau_;

    Eigen::Vector3d u_; //tau_ after projection

    void initializeMatrices();

    void forwardRecursion(Eigen::Vector3d q, Eigen::Vector3d dq, Eigen::Vector3d ddq);

    void backwardRecursion(Eigen::Vector3d q);



};

#endif
