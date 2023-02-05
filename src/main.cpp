
#include <vector>
#include <chrono>

#include <Eigen/Dense>

#include "dynamical_model.h"
#include "kinematic_model.h"
#include "inverse_kinematic.h"
#include "robot.h"


int main(int argc, char *argv[])
{

    const int DOFS = 6;
    //otherwise it returns always the same number
    srand((unsigned int) time(0));



    Eigen::VectorXd q_test;
    Eigen::VectorXd dq;
    Eigen::VectorXd ddq;

    q_test.resize(DOFS);
    dq.resize(DOFS);
    ddq.resize(DOFS);

    int first_link = 0;
    int last_link = 5;

    q_test <<  3.49677, 4.33261, 1.39203, 5.98386, 2.38322, 6.28319;
    //q_test.setZero();
    //q_test <<  0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    std::string param_file = "../robots/robot_parameters.json";
    Robot robot;
    robot.buildRobotFromFile(param_file);

    KinematicModel km(robot);
    km.setQ(q_test);
    km.computeJacobian();
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "--------------------Jacobian--------------------" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    Eigen::MatrixXd jacobian;
    km.getJacobian(jacobian);
    std::cout  << jacobian << std::endl;
    std::cout << "fk --> " << km.getTrans() << std::endl;
/*
    InverseKinematic ik;
    Eigen::Vector3d desired_pos; 
    desired_pos << -0.81725,-0.10915,-0.005491;
    Eigen::VectorXd solution;
    solution.resize(DOFS);
    ik.setDesiredPos(desired_pos);
    ik.setQk(q_test);
    ik.solveIk(solution);
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "--------------------ik solution--------------------" << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << solution << std::endl;
    km.setQ(solution);
    std::vector<int> link_origins = {first_link,last_link};
    km.computeForwardKinematic(link_origins);
    std::cout << "fk --> " << km.getTrans() << std::endl;
    std::cout << "Final cartisian error norm --> " << (desired_pos-km.getTrans()).norm() << std::endl;




    Eigen::VectorXd g;
    g.resize(DOFS);
    DynamicalModel model_g;
     Eigen::Vector3d gravity_g;
    gravity_g << 0,0,9.81;
    dq.setZero();
    ddq.setZero();
    auto start_g = std::chrono::high_resolution_clock::now();
    g = model_g.rnea(q_test, dq, ddq, gravity_g);
    auto end_g = std::chrono::high_resolution_clock::now();
    auto duration_g = std::chrono::duration_cast<std::chrono::microseconds>(end_g-start_g);

    DynamicalModel model_M;
    Eigen::MatrixXd M;
    M.resize(DOFS,DOFS);
    M.setZero();

    Eigen::Vector3d gravity;
    gravity << 0,0,0;

    Eigen::VectorXd M_i;
    M_i.resize(DOFS);
    M_i.setZero();
    dq.setZero();
    auto start_m = std::chrono::high_resolution_clock::now();
    for(short int i=0; i<DOFS; i++){
        ddq.setZero();
        ddq[i] = 1;
        M_i = model_M.rnea( q_test, dq, ddq, gravity);
        for(short int j=0; j<DOFS; j++){
            M(j,i) = M_i[j];
        }
        M_i.setZero();
    }
    auto end_m = std::chrono::high_resolution_clock::now();
    auto duration_m = std::chrono::duration_cast<std::chrono::microseconds>(end_m-start_m);


    std::cout << "--------------------g vector--------------------" << std::endl;
    std::cout << g << std::endl;
    std::cout << "Compute in " << duration_g.count() << std::endl;

    std::cout << "--------------------M matrix--------------------" << std::endl;
    std::cout << M << std::endl;
    std::cout << "Compute in " << duration_m.count() << std::endl;
*/
    return 0;

}
