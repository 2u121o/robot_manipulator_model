#ifndef DYNAMICAL_MODEL_H
#define DYNAMICAL_MODEL_H

#include <iostream>
#include <vector>
#include <math.h>

#include <Eigen/Dense>

#include "kinematic_model.h"
#include "robot.h"

class DynamicalModel{

public:
    DynamicalModel(Robot robot);
    ~DynamicalModel();

    Eigen::VectorXd rnea(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq, const Eigen::Vector3d gravity);

private:
    const int DOFS = 6;
    Eigen::VectorXd m_;

    KinematicModel kinematic_model_;
    Robot robot_;

    std::vector<Eigen::Matrix3d> inertia_; //3X3N tensor contains inertia matrix for each link
    std::vector<Eigen::Vector3d> cog_;     //3XN matrix contains cog vector for each link

    std::vector<Eigen::Vector3d> omega_;    //angular velocity
    std::vector<Eigen::Vector3d> d_omega_;  //angular acceleration
    std::vector<Eigen::Vector3d> a_;      //linear acceleration at frame i in frame i
    std::vector<Eigen::Vector3d> ac_;     //linear acceleration of the cog in frame i

    std::vector<Eigen::Vector3d> f_;        //linear forces applied at the origin of the RF
    std::vector<Eigen::Vector3d> tau_;      //tourques

    Eigen::VectorXd u_;                     //tau_ after projection

    void initializeMatrices(const Eigen::Vector3d gravity);

    void forwardRecursion(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq);

    void backwardRecursion(const Eigen::VectorXd &q);



};

#endif
