#ifndef KINEMATIC_MODEL_H
#define KINEMATIC_MODEL_H

#include <iostream>
#include <vector>
#include <forward_list>

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
template <typename T>
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
        KinematicModel(Robot<T> robot);

        /** @brief Compute the forward kinematics between two specified links 
        *  
        *  @param idx_links ...........
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
        void setQ(const Eigen::Matrix<T, Eigen::Dynamic, 1> &q);

        /** @brief Get the joint position q
        *
        *  @return is the actual joint position
        */
        Eigen::Matrix<T, Eigen::Dynamic, 1> getQ();

        /** @brief Get the position of the end effector in the base frame
        *
        *  @return a three dimensional vector that represent the positoin
        *          of the endeffector in the base frame
        */
        Eigen::Matrix<T, 3, 1> getTrans(const int start_link_idx, const int end_link_idx);

        /** @brief Get the orientation of the end effector w.r.t. the base frame
        *
        *  @return a 3X3 orthonormal matrix i.e., R^TR = RR^T = I follows R^T=R^{-1}, where 
        *                I is the identity matrix and T represents the transpose. This matrix
        *                describes the orientation of the endeffector with respect to the 
        *                base frame 
        */
        Eigen::Matrix<T, 3, 3> getR(const int start_link_idx, const int end_link_idx);

        /** @brief Get the geometric Jacobian from the base frame to the endeffector
        *
        *  @param jacobian is a 6xn matrix, where n is the number of the robots dofs. 
        *                  The matrix represent the geometric Jacobian from the robot
        *                  base to the endeffector
        */
        void getJacobian(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &jacobian);

    private:

    
        //! @brief Robots degree of freedom 
        int dofs_; 

        //! @brief Contains the kinematic and dynamic parameters of the robot
        Robot<T> robot_;

        //! @brief Joint position in radiants 
        Eigen::Matrix<T, Eigen::Dynamic,1> q_;

        //! @brief Three dimensional vector that represnts the translation between the origin of two reference frames
        Eigen::Matrix<T, Eigen::Dynamic,1> trans_;

        //! @brief Orthonormal matrix that represents the rotation between two reference frames
        Eigen::Matrix<T, Eigen::Dynamic,Eigen::Dynamic> R_;

        //! @brief Geometric jacobian between two referene frame. The dimension of the matrix is 6 X dofs_
        Eigen::Matrix<T, Eigen::Dynamic,Eigen::Dynamic> jacobian_;

        //! @brief Angle between the x axis of two reference frames about the z axis of the first reference frame.
        //! The angle is positive when the rotation is made conunter-clockwise.
        Eigen::Matrix<T, Eigen::Dynamic, 1> theta_;

        //! @brief Coordinate of two reference frame along the z axis of the previous reference frame
        Eigen::Matrix<T, Eigen::Dynamic, 1> d_;

        //! @brief Distance between the origin of two reference frame.
        Eigen::Matrix<T, Eigen::Dynamic, 1> a_;

        //! @brief Angle between the z axes of two reference frame about the x axis of the sencond 
        //! reference frame. The angle is positive when the rotation is made conunter-clockwise.
        Eigen::Matrix<T, Eigen::Dynamic, 1> alpha_;
        
        /**
         * @brief Resize all the necessary variable and initialize the kinematic  parameters with 
         * the information provided by the Robot object.
         * 
         */
        void initVariables();
        

};

template <typename T>
KinematicModel<T>::KinematicModel()
{

}

template <typename T>
KinematicModel<T>::KinematicModel(Robot<T> robot):robot_(robot)
{

    initVariables();
    // std::cout << "Kinematic model constructed! " << std::endl;
}

template <typename T>
void KinematicModel<T>::initVariables()
{

    std::forward_list<Link<T>> links;
    robot_.getLinks(links);

    dofs_=  robot_.getDofs();

    q_.resize(dofs_);
    jacobian_.resize(dofs_, dofs_);
    theta_.resize(dofs_);
    d_.resize(dofs_);
    a_.resize(dofs_);
    alpha_.resize(dofs_);

    q_.setZero();
    jacobian_.setZero();

    int i=0;  
    DHParams<T> dh_params;
    for(auto link : links)
    {
        link.getDHParams(dh_params);

        theta_(i) = dh_params.theta;
        d_(i) = dh_params.d;
        a_(i) = dh_params.a;
        alpha_(i) = dh_params.alpha;
        i++;

    }

    R_.resize(3*dofs_,3*dofs_);
    trans_.resize(3*dofs_);

    R_.setIdentity();
    trans_.setIdentity();

    //UR5 from https://www.universal-robots.com/articles/ur/application-installation/dh-parameters-for-calculations-of-kinematics-and-dynamics/
    //change also in computeForwardKinematic() if change robot
    // alpha_ << M_PI/2, 0, 0, M_PI/2, -M_PI/2, 0.0;
    // a_ << 0.0, -0.425, -0.39225, 0.0, 0.0, 0.0;
    // d_ << 0.089159, 0.0, 0.0, 0.10915, 0.09465, 0.0823;
}

template <typename T>
void KinematicModel<T>::computeForwardKinematic(const int start_link_idx, const int end_link_idx){
 

    Eigen::Matrix<T, 3, 3> R; 
    Eigen::Matrix<T, 3, 1> trans;

    R_.resize(3*dofs_,3*dofs_);
    trans_.resize(3*dofs_);

    R_.setIdentity();
    trans_.setIdentity();
   
    for(int i=start_link_idx; i<end_link_idx; ++i)
    { 
        double theta = theta_(i);
        double alpha = alpha_(i);
        double a = a_(i);
        double d = d_(i);
        
        theta += q_[i];

        if constexpr(std::is_same<T, double>::value)
        {
            double cos_theta, sin_theta;
            double cos_alpha, sin_alpha;

            sincos(theta, &sin_theta, &cos_theta);
            sincos(alpha, &sin_alpha, &cos_alpha);

            R << cos_theta, -sin_theta*cos_alpha, sin_theta*sin_alpha,
                 sin_theta, cos_theta*cos_alpha , -cos_alpha*sin_alpha,
                  0           , sin_alpha              , cos_alpha;
            trans << a*cos_theta, a*sin_theta, d;
        }
        else if constexpr(std::is_same<T, float>::value)
        {
            float cos_theta, sin_theta;
            float cos_alpha, sin_alpha;

            sincosf(theta, &sin_theta, &cos_theta);
            sincosf(alpha, &sin_alpha, &cos_alpha);

            R << cos_theta, -sin_theta*cos_alpha, sin_theta*sin_alpha,
                 sin_theta, cos_theta*cos_alpha , -cos_alpha*sin_alpha,
                  0           , sin_alpha              , cos_alpha;

            trans << a*cos_theta, a*sin_theta, d;
        }
        
        trans_.segment(i*3,+3) += R_.block(i*3,i*3,3,3)*trans ;
        R_.block(i*3,i*3,3,3) *= R;
        
    }

}

template <typename T>
void KinematicModel<T>::computeJacobian()
{

    Eigen::Matrix<T, 3, 1> pos_ee_absolute;
    pos_ee_absolute.setZero();

    computeForwardKinematic(0,dofs_-1);
    pos_ee_absolute = trans_;

    Eigen::Matrix<T, 3, 1> pos_ee_relative;
    pos_ee_relative.setZero();

    Eigen::Matrix<T, 3, 1> z_i_minus_one;
    z_i_minus_one.setZero();

    Eigen::Matrix<T, 3, 1> z_i_minus_one_const;
    z_i_minus_one_const << 0,0,1;

    for(int i=0; i<dofs_; i++){
        computeForwardKinematic(0,i);
        z_i_minus_one = R_*z_i_minus_one_const;
        pos_ee_relative = pos_ee_absolute - trans_;
        jacobian_.block(0,i,3,1) = z_i_minus_one.cross(pos_ee_relative);
        jacobian_.block(3,i,3,1) = z_i_minus_one;
    }
}

template <typename T>
void KinematicModel<T>::setQ(const Eigen::Matrix<T, Eigen::Dynamic, 1> &q)
{
    q_ = q;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> KinematicModel<T>::getQ()
{
    return q_;
}

template <typename T>
Eigen::Matrix<T, 3, 1> KinematicModel<T>::getTrans(const int start_link_idx, const int end_link_idx)
{
    return trans_.segment(start_link_idx*3,end_link_idx*3-start_link_idx*3);
}

template <typename T>
Eigen::Matrix<T, 3, 3> KinematicModel<T>::getR(const int start_link_idx, const int end_link_idx)
{
    
    return R_.block(start_link_idx*3,start_link_idx*3,end_link_idx*3-start_link_idx*3,end_link_idx*3-start_link_idx*3);
}

template <typename T>
void KinematicModel<T>::getJacobian(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &jacobian)
{
    jacobian = jacobian_;
}


#endif