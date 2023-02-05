#ifndef KINEMATIC_MODEL_H
#define KINEMATIC_MODEL_H

#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "robot.h"
#include "link.h"
#include "utils.h"

class KinematicModel{

    public:
        KinematicModel();
        KinematicModel(Robot robot);
        ~KinematicModel();

        /** @brief Compute the forward kinematics between two specified links 
        *  
        *  @param idx_links is the vector containing the indeces of the link
        *                   origins for which we want to compute the forward
        *                   kinematic
        */
        void computeForwardKinematic(const std::vector<int> idx_links);

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
        Eigen::Vector3d getTrans();

        /** @brief Get the orientation of the end effector w.r.t. the base frame
        *
        *  @return a 3X3 orthonormal matrix i.e., R^TR = RR^T = I follows R^T=R^{-1}, where 
        *                I is the identity matrix and T represents the transpose. This matrix
        *                describes the orientation of the endeffector with respect to the 
        *                base frame 
        */
        Eigen::Matrix3d getR();

        /** @brief Get the geometric Jacobian from the base frame to the endeffector
        *
        *  @param jacobian is a 6xn matrix, where n is the number of the robots dofs. 
        *                  The matrix represent the geometric Jacobian from the robot
        *                  base to the endeffector
        */
        void getJacobian(Eigen::MatrixXd &jacobian);

    private:

        int dofs_; 

        Robot robot_;

        Eigen::VectorXd q_;
        Eigen::Vector3d cartesina_pos_;
        Eigen::Vector3d trans_;

        Eigen::MatrixXd jacobian_;
        Eigen::Matrix3d R_;


        //parameters of the DH table

        Eigen::VectorXd theta_;
        Eigen::VectorXd d_;
        Eigen::VectorXd a_;
        Eigen::VectorXd alpha_;

        void initVariables();
        

};

#endif