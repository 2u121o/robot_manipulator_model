#ifndef DYNAMICAL_MODEL_H
#define DYNAMICAL_MODEL_H

#include <iostream>
#include <vector>
#include <math.h>

#include <Eigen/Dense>

#include "kinematic_model.h"
#include "robot.h"

/**
 * @class DynamicalModel
 * @brief  Dynamical model of the robot.
 *
 * This class porvides functionalities to comupte the dynamics of an 
 * arbitrary robot. 
 *
 * @startuml
 * class DynamicalModel {
 * +mass_: Eigen::VectorXd
 * +dofs_: int
 * +kinematic_model_: KinematicModel
 *}
 * @enddot
 * @relates KinematicModel
 */
template <typename T>
class DynamicalModel{

public:
    DynamicalModel(Robot<T> robot);

    Eigen::Matrix<T, Eigen::Dynamic, 1> rnea(const Eigen::Matrix<T, Eigen::Dynamic, 1> &q, const Eigen::Matrix<T, Eigen::Dynamic, 1> &dq, const Eigen::Matrix<T, Eigen::Dynamic, 1> &ddq, const Eigen::Matrix<T, 3, 1> gravity);

private:

    int dofs_;

    Eigen::Matrix<T, Eigen::Dynamic, 1> mass_;

    
    KinematicModel<T> kinematic_model_;
    Robot<T> robot_;

    //cogs expressed wrt the next frame so ready to use in the formulas without transformation
    std::vector<Eigen::Matrix<T, 3, 3>> inertia_; //3X3N tensor contains inertia matrix for each link
    std::vector<Eigen::Matrix<T, 3, 1>> cog_;     //3XN matrix contains cog vector for each link

    std::vector<Eigen::Matrix<T, 3, 1>> omega_;    //angular velocity
    std::vector<Eigen::Matrix<T, 3, 1>> d_omega_;  //angular acceleration
    std::vector<Eigen::Matrix<T, 3, 1>> a_;      //linear acceleration at frame i in frame i
    std::vector<Eigen::Matrix<T, 3, 1>> ac_;     //linear acceleration of the cog in frame i

    std::vector<Eigen::Matrix<T, 3, 1>> f_;        //linear forces applied at the origin of the RF
    std::vector<Eigen::Matrix<T, 3, 1>> tau_;      //tourques

    Eigen::Matrix<T, Eigen::Dynamic, 1> u_;                     //tau_ after projection

    void initializeMatrices(const Eigen::Matrix<T, 3, 1> gravity);

    void forwardRecursion(const Eigen::Matrix<T, Eigen::Dynamic, 1> &q, const Eigen::Matrix<T, Eigen::Dynamic, 1> &dq, const Eigen::Matrix<T, Eigen::Dynamic, 1> &ddq);

    void backwardRecursion(const Eigen::Matrix<T, Eigen::Dynamic, 1> &q);

    /**
     * @brief Resize vectors using the dofs in the robot object.
     * 
     */
    void resizeVariables();

    /**
     * @brief Initializes e the vectors containing the dynaimc parameters with the value from robot object
     * 
     */
    void initializeDynamicParameters();

};


template <typename T>
DynamicalModel<T>::DynamicalModel(Robot<T> robot):robot_(robot)
{
    dofs_ = robot.getDofs();
    
    resizeVariables();
    initializeDynamicParameters();
    kinematic_model_ = KinematicModel<T>(robot_);
}


template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> DynamicalModel<T>::rnea(const Eigen::Matrix<T, Eigen::Dynamic, 1> &q, const Eigen::Matrix<T, Eigen::Dynamic, 1> &dq, const Eigen::Matrix<T, Eigen::Dynamic, 1> &ddq, const Eigen::Matrix<T, 3, 1> gravity)
{
    initializeMatrices(gravity);

    kinematic_model_.setQ(q);
    kinematic_model_.computeForwardKinematic(0,q.size());

    forwardRecursion(q, dq, ddq);
    backwardRecursion(q);

    return  u_;
}

template <typename T>
void DynamicalModel<T>::forwardRecursion(const Eigen::Matrix<T, Eigen::Dynamic, 1> &q, const Eigen::Matrix<T, Eigen::Dynamic, 1> &dq, const Eigen::Matrix<T, Eigen::Dynamic, 1> &ddq)
{
    Eigen::Matrix<T, 3, 1> z;
    z << 0, 0, 1;

    Eigen::Matrix<T, 3, 3> R;
    Eigen::Matrix<T, 3, 3> R_transpose;
    Eigen::Matrix<T, 3, 1> t;
  
    for(short int i=0; i<dofs_; i++){

        R = kinematic_model_.getR(i,i+1);
        t = kinematic_model_.getTrans(i,i+1);

        R_transpose = R.transpose();

        omega_.at(i+1) = R_transpose*(omega_.at(i) + dq[i]*z);
        d_omega_.at(i+1) = R_transpose*(d_omega_.at(i) + ddq[i]*z + dq[i]*(omega_.at(i).cross(z)));
        a_.at(i+1) = R_transpose*a_.at(i) + d_omega_.at(i+1).cross(R_transpose*t) + omega_.at(i+1).cross(omega_.at(i+1).cross(R_transpose*t));    
        ac_.at(i) = a_.at(i+1) + d_omega_.at(i+1).cross(cog_.at(i)) + omega_.at(i+1).cross(omega_.at(i+1).cross(cog_.at(i)));
    }

}

template <typename T>
void DynamicalModel<T>::backwardRecursion(const Eigen::Matrix<T, Eigen::Dynamic, 1> &q)
{

    Eigen::Matrix<T, 3, 1> z;
    z << 0, 0, 1;

    Eigen::Matrix<T, 3, 3> R;
    Eigen::Matrix<T, 3, 3> R_transpose;
    Eigen::Matrix<T, 3, 1> t;
    Eigen::Matrix<T, 3, 3> Rp1;

    R.setZero();
    t.setZero();

    //the first transfomation is identity because it transfom
    //from the tip to the environment so it remain like it is
    Rp1 << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;

    for(short int i=dofs_-1; i>=0; i--){

        // kinematic_model_.setQ(q);
        R = kinematic_model_.getR(i,i+1);
        t = kinematic_model_.getTrans(i,i+1);

        R_transpose = R.transpose();

        f_.at(i) = Rp1*f_.at(i+1) + mass_[i]*(ac_.at(i));
       
        tau_.at(i) = Rp1*tau_.at(i+1) + (Rp1*f_.at(i+1)).cross(cog_.at(i)) - f_.at(i).cross(R_transpose*t + cog_.at(i)) +
                    inertia_.at(i)*(d_omega_.at(i+1)) + omega_.at(i+1).cross(inertia_.at(i) * (omega_.at(i+1)));

        u_[i] = tau_.at(i).transpose()*R_transpose*z;
        Rp1 = R;
    }

}

template <typename T>
void DynamicalModel<T>::initializeMatrices(const Eigen::Matrix<T, 3, 1> gravity)
{
    for(short int i=0; i<dofs_+1; i++){
        omega_.at(i).setZero();
        d_omega_.at(i).setZero();
        a_.at(i).setZero();
        f_.at(i).setZero();
        tau_.at(i).setZero();
        if(i<dofs_) ac_.at(i).setZero();
    }
    a_.at(0) = -gravity;

    u_.setZero();   
}

template <typename T>
void DynamicalModel<T>::resizeVariables()
{
    mass_.resize(dofs_);
    inertia_.resize(dofs_);
    cog_.resize(dofs_);

    omega_.resize(dofs_+1);
    d_omega_.resize(dofs_+1);
    a_.resize(dofs_+1);
    ac_.resize(dofs_);

    f_.resize(dofs_+1);
    tau_.resize(dofs_+1);

    u_.resize(dofs_);
}

template <typename T>
void DynamicalModel<T>::initializeDynamicParameters()
{
    std::forward_list<Link<T>> links;
    robot_.getLinks(links);
    int link_id=0;
    DynamicParameters<T> dynamic_parameters;

    for(auto link : links)
    {
        link.getDynamicParameters(dynamic_parameters);

        mass_(link_id) = dynamic_parameters.mass;
        cog_.at(link_id) = dynamic_parameters.com;

        Eigen::Matrix<T, 6, 6> inertia = dynamic_parameters.inertia;
    
        int count = 0;
        const int MATRIX_DIM = 3;
        for(int i=0; i<MATRIX_DIM; i++){
            for(int j=i; j<MATRIX_DIM; j++){
                double element = inertia(count++);
                inertia_.at(i)(i,j) = element;
                inertia_.at(i)(j,i) = element;
            }
        }
        link_id++;
    }
}

#endif
