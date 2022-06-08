
#include "dynamical_model.h"

#include <Eigen/Dense>

int main(int argc, char *argv[])
{
    //otherwise it returns always the same number
    srand((unsigned int) time(0));

    DynamicalModel model_n;
    Eigen::VectorXd q;
    Eigen::VectorXd dq;
    Eigen::VectorXd ddq;

    q.resize(3);
    dq.resize(3);
    ddq.resize(3);

    q.setRandom();


    Eigen::VectorXd n;
    n.resize(3);
    dq.setRandom();
    ddq.setZero();
    n = model_n.rnea(q, dq, ddq);

    DynamicalModel model_M;
    Eigen::MatrixXd M;
    M.resize(3,3);
    M.setZero();

    Eigen::VectorXd M_i;
    M_i.resize(3);
    M_i.setZero();
    dq.setZero();
    for(short int i=0; i<3; i++){
        ddq.setZero();
        ddq[i] = 1;
        M_i = model_M.rnea( q, dq, ddq);
        for(short int j=0; j<3; j++){
            M(j,i) = M_i[j];
        }
        M_i.setZero();
    }


    std::cout << "--------------------n vector--------------------" << std::endl;
    std::cout << n << std::endl;

    std::cout << "--------------------M matrix--------------------" << std::endl;
    std::cout << M << std::endl;

    return 0;

}
