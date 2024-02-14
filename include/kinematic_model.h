#ifndef KINEMATIC_MODEL_H
#define KINEMATIC_MODEL_H

#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

#include "robot.h"
#include "link.h"
#include "utils.h"

/**
 * @class KinematicModel
 * @brief  Kinematic model of the robot.
 *
 * This class porvides functionalities to comupte the geometric Jacobian between two 
 * arbitrary frames and also the forward kinematic. 
 *
 * @startuml
 * class KinematicModel {
 * +dofs_: int
 * +cognome: string
 *}
 * @enddot
 */
class KinematicModel{

    public:


        /**
         * @brief Construct a new Kinematic Model object
         * 
         */
        KinematicModel();

        /**
         * @brief Construct a new Kinematic Model object with the information of 
         * robot contained in the Robot class
         * 
         * @param robot contains the kinematics and dynamics parameters of the 
         * robot
         */
        KinematicModel(Robot robot);

        /**
         * @brief Destroy the Kinematic Model object
         * 
         */
        ~KinematicModel();

        /** @brief Compute the forward kinematics between two specified links 
        *  
        *  @param idx_links is the vector containing the indeces of the link
        *                   origins for which we want to compute the forward
        *                   kinematic
        */
        void computeForwardKinematic(const int start_link_idx, const int end_link_idx);

        /** @brief Compute the Jacobian at the endeffector starting from the 
        *         robot base
        */
        void computeJacobian();

        /** @brief Set the joint position q
        *
        *  @param q is the joint position
        */
        void setQ(const Eigen::VectorXd &q);

        /** @brief Get the joint position q
        *
        *  @return is the actual joint position
        */
        Eigen::VectorXd getQ();

        /** @brief Get the position of the end effector in the base frame
        *
        *  @return a three dimensional vector that represent the positoin
        *          of the endeffector in the base frame
        */
        Eigen::Vector3d getTrans(const int start_link_idx, const int end_link_idx);

        /** @brief Get the orientation of the end effector w.r.t. the base frame
        *
        *  @return a 3X3 orthonormal matrix i.e., R^TR = RR^T = I follows R^T=R^{-1}, where 
        *                I is the identity matrix and T represents the transpose. This matrix
        *                describes the orientation of the endeffector with respect to the 
        *                base frame 
        */
        Eigen::Matrix3d getR(const int start_link_idx, const int end_link_idx);

        /** @brief Get the geometric Jacobian from the base frame to the endeffector
        *
        *  @param jacobian is a 6xn matrix, where n is the number of the robots dofs. 
        *                  The matrix represent the geometric Jacobian from the robot
        *                  base to the endeffector
        */
        void getJacobian(Eigen::MatrixXd &jacobian);

    private:
    
    
        //! @brief Robots degree of freedom 
        int dofs_; 

        //! @brief Contains the kinematic and dynamic parameters of the robot
        Robot robot_;

        //! @brief Joint position in radiants 
        Eigen::VectorXd q_;

        //! @brief Three dimensional vector that represnts the translation between the origin of two reference frames
        Eigen::VectorXd trans_;

        //! @brief Orthonormal matrix that represents the rotation between two reference frames
        Eigen::MatrixXd R_;

        //! @brief Geometric jacobian between two referene frame. The dimension of the matrix is 6 X dofs_
        Eigen::MatrixXd jacobian_;

        //! @brief Angle between the x axis of two reference frames about the z axis of the first reference frame.
        //! The angle is positive when the rotation is made conunter-clockwise.
        Eigen::VectorXd theta_;

        //! @brief Coordinate of two reference frame along the z axis of the previous reference frame
        Eigen::VectorXd d_;

        //! @brief Distance between the origin of two reference frame.
        Eigen::VectorXd a_;

        //! @brief Angle between the z axes of two reference frame about the x axis of the sencond 
        //! reference frame. The angle is positive when the rotation is made conunter-clockwise.
        Eigen::VectorXd alpha_;

        float cos_theta;
        float sin_theta;
        float cos_alpha;
        float sin_alpha;

        float theta_prev_;
        float alpha_prev_;
        
        /**
         * @brief Resize all the necessary variable and initialize the kinematic  parameters with 
         * the information provided by the Robot object.
         * 
         */
        void initVariables();
        

};

#endif