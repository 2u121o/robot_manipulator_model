#ifndef DYNAMICAL_MODEL_H
#define DYNAMICAL_MODEL_H

#include <vector>

#include <Eigen/Dense>

class DynamicalModel{

public:
    DynamicalModel();

    void forwardRecursion(Eigen::VectorXd q, Eigen::VectorXd dq, Eigen::VectorXd ddq);

private:
    const int DOFS = 3;
    Eigen::VectorXd m_;

    std::vector<Eigen::Matrix3d> inertia_; //3X3N tensor contains inertia matrix for each link
    std::vector<Eigen::Vector3d> cog_;     //3XN matrix contains cog vector for each link

};

#endif
