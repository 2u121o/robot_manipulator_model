
#include <vector>
#include <chrono>

#include <Eigen/Dense>

#include "dynamical_model.h"
#include "kinematic_model.h"
#include "inverse_kinematic.h"
#include "robot.h"


int main(int argc, char *argv[])
{
    // TEST
    // sudo perf stat  ./robot_maipulator_model
    // sudo perf stat -B -e cache-references,cache-misses,cycles,instructions,branches,faults,migrations ./robot_maipulator_model

    // sudo perf record -e cache-misses ./robot_maipulator_model
    // sudo perf report -v

    // valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all  ./robot_maipulator_model


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

    //q_test <<  3.49677, 4.33261, 1.39203, 5.98386, 2.38322, 6.28319;
    q_test.setZero();
    //q_test <<  0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    std::string param_file = "../robots/robot_parameters.json";
    Robot robot;
    robot.buildRobotFromFile(param_file);

    // std::forward_list<Link> links;
    // robot.getLinks(links);

    // DynamicParameters dynamic_parameter;
    // for(Link& link : links)
    // {
    //     link.getDynamicParameters(dynamic_parameter);
    //     std::cout << "dynamic_parameter.mass: " <<  dynamic_parameter.mass << std::endl;
    // }

    Eigen::VectorXd g;
    g.resize(DOFS);
    DynamicalModel model_g(robot);
    Eigen::Vector3d gravity_g;
    gravity_g << 0,0, -9.81;
    dq.setZero();
    ddq.setZero();
    q_test.setZero();
    q_test(1) = 1.0;
    q_test(2) = 1.0;
    q_test(3) = 1.0;

    int NUM_ITERATIONS = 1000000;
    auto duration_mean_g  = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::duration<double>(0.0));
    double duration_max = 0.0;

    for(int i=0; i<NUM_ITERATIONS; i++)
    {
        auto start_g = std::chrono::high_resolution_clock::now();
        g = model_g.rnea(q_test, dq, ddq, gravity_g);
        auto end_g = std::chrono::high_resolution_clock::now();
        auto duration_g = std::chrono::duration_cast<std::chrono::microseconds>(end_g-start_g);
        duration_mean_g += duration_g;
        if(duration_max<duration_g.count())
        {
            duration_max = duration_g.count();
        }
    }
    double mean_computation_g = duration_mean_g.count()/static_cast<double>(NUM_ITERATIONS);

    std::cout << "--------------------g vector--------------------" << std::endl;
    std::cout << g << std::endl;
    std::cout << "Compute in (mean) " << mean_computation_g << " us" << std::endl;
    std::cout << "Max computation time " << duration_max << " us" << std::endl;

    //EXPECTED
     // --------------------g vector--------------------
    // -4.66294e-15
    //     -156.408
    //     -50.976
    //     -1.10659
    //     -1.65226
    // -0.223362
    // Compute in 19

    // KinematicModel km(robot);
    // km.setQ(q_test);
    // km.computeJacobian();
    // std::cout << "------------------------------------------------" << std::endl;
    // std::cout << "--------------------Jacobian--------------------" << std::endl;
    // std::cout << "------------------------------------------------" << std::endl;
    // Eigen::MatrixXd jacobian;
    // km.getJacobian(jacobian);
    // std::cout  << jacobian << std::endl;
    // std::cout << "fk --> " << km.getTrans() << std::endl;

    // InverseKinematic ik(robot);
    // Eigen::Vector3d desired_pos; 
    // desired_pos << -0.41725,-0.10915,-0.005491;
    // Eigen::VectorXd solution;
    // solution.resize(DOFS);
    // ik.setDesiredPos(desired_pos);
    // ik.setQk(q_test);
    // ik.solveIk(solution);
    // std::cout << "---------------------------------------------------" << std::endl;
    // std::cout << "--------------------ik solution--------------------" << std::endl;
    // std::cout << "---------------------------------------------------" << std::endl;
    // std::cout << solution << std::endl;
    // km.setQ(solution);
    // std::vector<int> link_origins = {first_link,last_link};
    // km.computeForwardKinematic(link_origins);
    // std::cout << "fk --> " << km.getTrans() << std::endl;
    // std::cout << "Final cartisian error norm --> " << (desired_pos-km.getTrans()).norm() << std::endl;




    // DynamicalModel model_M(robot);
    // Eigen::MatrixXd M;
    // M.resize(DOFS,DOFS);
    // M.setZero();

    // Eigen::Vector3d gravity;
    // gravity << 0,0,0;

    // Eigen::VectorXd M_i;
    // M_i.resize(DOFS);
    // M_i.setZero();
    // dq.setZero();
    // auto start_m = std::chrono::high_resolution_clock::now();
    // for(short int i=0; i<DOFS; i++){
    //     ddq.setZero();
    //     ddq[i] = 1;
    //     M_i = model_M.rnea( q_test, dq, ddq, gravity);
    //     for(short int j=0; j<DOFS; j++){
    //         M(j,i) = M_i[j];
    //     }
    //     M_i.setZero();
    // }
    // auto end_m = std::chrono::high_resolution_clock::now();
    // auto duration_m = std::chrono::duration_cast<std::chrono::microseconds>(end_m-start_m);


   
    // std::cout << "--------------------M matrix--------------------" << std::endl;
    // std::cout << M << std::endl;
    // std::cout << "Compute in " << duration_m.count() << std::endl;

    // q_test.setZero();
    // q_test(1) = 1.0;
    // q_test(2) = 1.0;
    // q_test(3) = 1.0;
    // --------------------g vector--------------------
    // -4.66294e-15
    //     -156.408
    //     -50.976
    //     -1.10659
    //     -1.65226
    // -0.223362
    // Compute in 19
    // --------------------M matrix--------------------
    // 13.7933  -0.232461 -0.0234105  -0.192564   0.346809  0.0724387
    // -0.232533    14.4179    5.65443   0.107887   0.237748  0.0338296
    // -0.0234822    5.65443    3.69536  0.0677056   0.161394  0.0256859
    // -0.192464   0.107715  0.0675336  0.0170558  0.0154857 0.00439645
    // 0.346809   0.237748   0.161394  0.0154857  0.0809201  0.0185367
    // 0.0724387  0.0338296  0.0256859 0.00439645  0.0185367  0.0044264
    // Compute in 76

    return 0;

}
