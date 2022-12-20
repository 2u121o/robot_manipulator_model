
#include "dynamical_model.h"

#include <Eigen/Dense>

int main(int argc, char *argv[])
{
    //otherwise it returns always the same number
    srand((unsigned int) time(0));

    DynamicalModel model_n;
    Eigen::Vector3d q;
    Eigen::Vector3d dq;
    Eigen::Vector3d ddq;

    //q.setRandom();
    q << 0, 0, M_PI/2;

    Eigen::Vector3d n;
    dq.setRandom();
    ddq.setZero();
    n = model_n.rnea(q, dq, ddq);

    Eigen::Vector3d g;
    dq.setZero();
    ddq.setZero();
    g = model_n.rnea(q, dq, ddq);

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


    std::cout << "--------------------n vector--------------------" << std::endl;
    std::cout << n << std::endl;

    std::cout << "--------------------g vector--------------------" << std::endl;
    std::cout << g << std::endl;

    std::cout << "--------------------M matrix--------------------" << std::endl;
    std::cout << M << std::endl;


    double I1zz, I2yy, I3xx, I3yy, I3zz;
    I1zz = 0.1;
    I2yy = 0.1;
    I3xx = 0.1;
    I3yy = 0.1;
    I3zz = 0.1;

    double m1, m2, m3;
    m1 = 1.0;
    m2 = 1.0;
    m3 = 1.0;

    double d1, d3;
    d1 = 0.2;
    d3 = 0.2;


    double l1 = 0.6;

    double m11 = I1zz+m1*std::pow(d1,2)+I2yy+(m2+m3)*std::pow(l1,2)+m3*std::pow(d3,2)*std::pow(cos(q(2)),2)
                +2*m3*l1*d3*cos(q(1))*cos(q(2))+(I3xx*std::pow(sin(q(2)),2)+I3yy*std::pow(cos(q(2)),2));

    double m12 = I2yy+m3*std::pow(d3,2)*std::pow(cos(q(2)),2)+m3*l1*d3*cos(q(1))*cos(q(2))+
                (I3xx*std::pow(sin(q(2)),2)+I3yy*std::pow(cos(q(2)),2));

    double m13 = -m3*l1*d3*sin(q(1))*sin(q(2));

    double m22 = I2yy+m3*std::pow(d3,2)*std::pow(cos(q(2)),2)+(I3xx*std::pow(sin(q(2)),2)+I3yy*std::pow(cos(q(2)),2));

    double m33 = I3zz+m3*std::pow(d3,2);
    
    Eigen::MatrixXd M_lag;
    M_lag.resize(3,3);
    M_lag.setZero();

    M_lag(0,0) = m11;
    M_lag(0,1) = m12;
    M_lag(0,2) = m13;
    M_lag(1,0) = m12;
    M_lag(1,1) = m22;
    M_lag(2,0) = m13;
    M_lag(2,2) = m33;

    std::cout << "--------------------M_lag matrix--------------------" << std::endl;
    std::cout << M_lag << std::endl;

    return 0;

}
