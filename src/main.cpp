//IMPLEMENTA UN 2R ROBOT CHE E MOLTO PIU SEMPLICE DA DEBBUGGARE
#include "dynamical_model.h"

#include <Eigen/Dense>

int main(int argc, char *argv[])
{
    //I set the seams for the random number generator
    //otherwise it returns always the same number
    srand((unsigned int) time(0));

   // DynamicalModel model_n;

    Eigen::Vector3d q;
    Eigen::Vector3d dq;
    Eigen::Vector3d ddq;

    q.setRandom();


    Eigen::Vector3d n;
    dq.setRandom();
    ddq.setZero();
   // model_n.rnea(n, q, dq, ddq);

    DynamicalModel model_M;
    Eigen::MatrixXd M;
    M.resize(3,3);
    M.setZero();

    Eigen::Vector3d M_i;
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


    //there is an error because the inertia matrix is not symmetric
    std::cout << "--------------------M matrix--------------------" << std::endl;
    std::cout << M << std::endl;

    return 0;

}
