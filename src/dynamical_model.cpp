#include "dynamical_model.h"

DynamicalModel::DynamicalModel(Robot robot):robot_(robot){

    dofs_ = robot.getDofs();
    
    resizeVariables();
    initializeDynamicParameters();
    kinematic_model_ = KinematicModel(robot_);
}



Eigen::VectorXd DynamicalModel::rnea(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq, const Eigen::Vector3d gravity){

    initializeMatrices(gravity);

    kinematic_model_.setQ(q);
    kinematic_model_.computeForwardKinematic(0,q.size());

    forwardRecursion(q, dq, ddq);
    backwardRecursion(q);

    return  u_;

}

void DynamicalModel::forwardRecursion(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq){

    Eigen::Vector3d z;
    z << 0, 0, 1;

    Eigen::Matrix3d R;
    Eigen::Matrix3d R_transpose;
    Eigen::Vector3d t;
  
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


void DynamicalModel::backwardRecursion(const Eigen::VectorXd &q){


    Eigen::Vector3d z;
    z << 0, 0, 1;

    Eigen::Matrix3d R;
    Eigen::Matrix3d R_transpose;
    Eigen::Vector3d t;
    Eigen::Matrix3d Rp1;

    R.setZero();
    t.setZero();

    //the first transfomation is identity because it transfom
    //from the tip to the environment so it remain like it is
    Rp1 << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;

    for(short int i=dofs_-1; i>=0; i--){

        kinematic_model_.setQ(q);
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

void DynamicalModel::initializeMatrices(const Eigen::Vector3d gravity){

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

void DynamicalModel::resizeVariables(){

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

void DynamicalModel::initializeDynamicParameters(){

    std::vector<Link> links;
    robot_.getLinks(links);
    for(int link_id=0; link_id<dofs_; link_id++){
         
        DynamicParameters dynamic_parameters;
        links[link_id].getDynamicParameters(dynamic_parameters);

        mass_(link_id) = dynamic_parameters.mass;
        cog_.at(link_id) = dynamic_parameters.com;

        Eigen::VectorXd iv = dynamic_parameters.inertia;
    
        int count = 0;
        const int MATRIX_DIM = 3;
        for(int i=0; i<MATRIX_DIM; i++){
            for(int j=i; j<MATRIX_DIM; j++){
                double element = iv(count++);
                inertia_.at(i)(i,j) = element;
                inertia_.at(i)(j,i) = element;
            }
        }

    }
}

DynamicalModel::~DynamicalModel(){}


