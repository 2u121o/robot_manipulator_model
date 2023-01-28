#ifndef KINEMATIC_MODEL_H
#define KINEMATIC_MODEL_H

#include <iostream>
#include <Eigen/Dense>

class KinematicModel{

    public:
        KinematicModel();
        ~KinematicModel();

        //TODO: 
        //      computeInverseKinematic
        //      computeJacobian

        void computeForwardKinematic(int idx_first_link, int idx_last_link);
        void computeJacobian();

        //Setter--------------------------
        void setQ(const Eigen::VectorXd &q);

        //Getter--------------------------
        Eigen::VectorXd getQ();
        Eigen::Vector3d getTrans();
        Eigen::Matrix3d getR();
        Eigen::MatrixXd getJacobian();

    private:

        const int DOFS = 6;

        Eigen::VectorXd q_;
        Eigen::Vector3d cartesina_pos_;
        Eigen::Vector3d trans_;

        Eigen::MatrixXd jacobian_;
        Eigen::Matrix3d R_;


        //parameters of the DH table
        Eigen::VectorXd alpha_;
        Eigen::VectorXd a_;
        Eigen::VectorXd d_;

        void initVariables();
        

};

#endif