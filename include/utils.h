#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>

/**
 * @brief Structure that contains the robots DH parameters
 * 
 */
template <typename T>
struct DHParams
{

    /**
     * @brief Angle between the x axis of two reference frames about the z axis of the first reference frame.
     * The angle is positive when the rotation is made counter-clockwise.
     * 
     */
    T theta;

    /**
     * @brief Coordinate of two reference frame along the z axis of the previous reference frame
     * 
     */
    T d;

    /**
     * @brief Distance between the origin of two reference frame.
     * 
     */
    T a;

    /**
     * @brief Angle between the z axes of two reference frame about the x axis of the second 
     * reference frame. The angle is positive when the rotation is made counter-clockwise.
     * 
     */
    T alpha;
};

/**
 * @brief Structure that contains the robots dynamic parameters.
 * 
 */
template <typename T>
struct DynamicParameters
{

    /**
     * @brief Links mass.
     * 
     */
    T mass;

    /**
     * @brief 3D vector that represents the position of the center of mass
     * in the reference frame of the current link.
     * 
     */
    Eigen::Matrix<T, 3, 1> com;

    /**
     * @brief 6D vector containing the component of the inertia tensor. The remaining
     * three components can be deduce since the matrix is symmetric.
     * 
     */
    Eigen::Matrix<T, 6, 1> inertia;
};

#endif