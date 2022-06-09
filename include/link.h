#ifndef LINK_H
#define LINK_H

#include <iostream>

#include <Eigen/Dense>

class Link{

public:
    Link(std::string name, double mass, Eigen::VectorXd inertia, Eigen::Vector3d cog, Eigen::Vector3d dh);
    Link(std::string name);

    void setName(std::string name);
    void setMass(double mass);
    void setInertia(Eigen::VectorXd inertia);
    void setCog(Eigen::Vector3d cog);
    void setDh(Eigen::Vector3d dh);

    std::string getName();
    double getMass();
    Eigen::VectorXd getInertia();
    Eigen::Vector3d getCog();
    Eigen::Vector3d getDh();

private:
    std::string name_;
    double mass_;
    Eigen::VectorXd inertia_;
    Eigen::Vector3d cog_;
    Eigen::Vector3d dh_;          //dh params, just d a alpha, not theta


};

#endif // LINK_H
