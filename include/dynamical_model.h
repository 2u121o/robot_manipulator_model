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
class DynamicalModel{

public:
    DynamicalModel(Robot robot);
    ~DynamicalModel();

    Eigen::VectorXd rnea(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq, const Eigen::Vector3d gravity);

private:

    int dofs_;

    Eigen::VectorXd mass_;

    
    KinematicModel kinematic_model_;
    Robot robot_;

    //cogs expressed wrt the next frame so ready to use in the formulas without transformation
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

#endif
