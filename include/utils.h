#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>

struct DHParams{
    double theta;
    double d;
    double a;
    double alpha;
};

struct DynamicParameters{
    double mass;
    Eigen::Vector3d com;
    Eigen::VectorXd inertia;
};

#endif