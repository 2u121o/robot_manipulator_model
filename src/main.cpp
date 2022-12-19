
#include "dynamical_model.h"

#include <Eigen/Dense>

#include <chrono>

int main(int argc, char *argv[])
{
    //otherwise it returns always the same number
    srand((unsigned int) time(0));

    
    Eigen::Vector3d q;
    Eigen::Vector3d dq;
    Eigen::Vector3d ddq;
    
    q.setRandom();


    Eigen::Vector3d n;
    dq.setRandom();
    ddq.setZero();

    DynamicalModel model_n;
    auto start_n = std::chrono::high_resolution_clock::now();
    n = model_n.rnea(q, dq, ddq);
    auto stop_n = std::chrono::high_resolution_clock::now();
    auto duration_n = std::chrono::duration_cast<std::chrono::microseconds>(stop_n-start_n);


    DynamicalModel model_M;
    Eigen::MatrixXd M;
    M.resize(3,3);
    M.setZero();

    Eigen::Vector3d M_i;
    M_i.setZero();
    dq.setZero();
     auto start_M = std::chrono::high_resolution_clock::now();
    for(short int i=0; i<3; i++){
        ddq.setZero();
        ddq[i] = 1;
        M_i = model_M.rnea( q, dq, ddq);
        for(short int j=0; j<3; j++){
            M(j,i) = M_i[j];
        }
        M_i.setZero();
    }
    auto stop_M = std::chrono::high_resolution_clock::now();
    auto duration_M = std::chrono::duration_cast<std::chrono::microseconds>(stop_M-start_M);

    std::cout << "--------------------n vector--------------------" << std::endl;
    std::cout << n << std::endl;
   std::cout << "computed in " << duration_n.count() << " microseconds" << std::endl;

    std::cout << "--------------------M matrix--------------------" << std::endl;
    std::cout << M << std::endl;
    std::cout << "computed in " << duration_M.count() << " microseconds" << std::endl;

    

return 0;
}
