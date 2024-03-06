#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>

/**
 * @brief Structure that contains the robots DH parameters
 * 
 */
struct DHParams{

    /**
     * @brief Angle between the x axis of two reference frames about the z axis of the first reference frame.
     * The angle is positive when the rotation is made counter-clockwise.
     * 
     */
    double theta;

    /**
     * @brief Coordinate of two reference frame along the z axis of the previous reference frame
     * 
     */
    double d;

    /**
     * @brief Distance between the origin of two reference frame.
     * 
     */
    double a;

    /**
     * @brief Angle between the z axes of two reference frame about the x axis of the second 
     * reference frame. The angle is positive when the rotation is made counter-clockwise.
     * 
     */
    double alpha;
};

/**
 * @brief Structure that contains the robots dynamic parameters.
 * 
 */
struct DynamicParameters{

    /**
     * @brief Links mass.
     * 
     */
    double mass;

    /**
     * @brief 3D vector that represents the position of the center of mass
     * in the reference frame of the current link.
     * 
     */
    Eigen::Vector3d com;

    /**
     * @brief 6D vector containing the component of the inertia tensor. The remaining
     * three components can be deduce since the matrix is symmetric.
     * 
     */
    Eigen::VectorXd inertia;
};

#endif