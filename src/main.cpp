
#include "dynamical_model.h"

#include <Eigen/Dense>

#include <chrono>

int main(int argc, char *argv[])
{

    const int DOFS = 6;
    //otherwise it returns always the same number
    srand((unsigned int) time(0));

    // DynamicalModel model_n;
    Eigen::VectorXd q;
    Eigen::VectorXd dq;
    Eigen::VectorXd ddq;

    q.resize(DOFS);
    dq.resize(DOFS);
    ddq.resize(DOFS);

    q <<  0.1, 0.5, 0.5, 0.3, 0.2, 0.3;
    //q.setRandom();
   // q << 0, 0, 0.0, 0.0, 0.0, 0;
    //q << 0, 0, 0, 0, 0, 0;
    //q << 0, 0.50, 0.50;


    // Eigen::VectorXd n;
    // n.resize(DOFS);
    // //dq.setRandom();
    // dq << 0,0,0,0,0,0;
    // ddq.setZero();
    // auto start = std::chrono::high_resolution_clock::now();
    // n = model_n.rnea(q, dq, ddq);
    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start);

    Eigen::VectorXd g;
    g.resize(DOFS);
    DynamicalModel model_g;
     Eigen::Vector3d gravity_g;
    gravity_g << 0,0,9.81;
    dq.setZero();
    ddq.setZero();
    auto start_g = std::chrono::high_resolution_clock::now();
    g = model_g.rnea(q, dq, ddq, gravity_g);
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
        M_i = model_M.rnea( q, dq, ddq, gravity);
        for(short int j=0; j<DOFS; j++){
            M(j,i) = M_i[j];
        }
        M_i.setZero();
    }
    auto end_m = std::chrono::high_resolution_clock::now();
    auto duration_m = std::chrono::duration_cast<std::chrono::microseconds>(end_m-start_m);

    // std::cout << "--------------------n vector--------------------" << std::endl;
    // std::cout << n << std::endl;
    // std::cout << "Compute in " << duration.count() << std::endl;

    std::cout << "--------------------g vector--------------------" << std::endl;
    std::cout << g << std::endl;
    std::cout << "Compute in " << duration_g.count() << std::endl;

    std::cout << "--------------------M matrix--------------------" << std::endl;
    std::cout << M << std::endl;
    std::cout << "Compute in " << duration_m.count() << std::endl;

    return 0;

}
