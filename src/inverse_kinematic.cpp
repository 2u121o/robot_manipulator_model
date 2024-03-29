#include "inverse_kinematic.h"

InverseKinematic::InverseKinematic(Robot robot):robot_(robot){
    initVariables();
}

void InverseKinematic::solveIk(Eigen::VectorXd &solution){

    kinematic_model_ = KinematicModel(robot_);
    kinematic_model_.setQ(q_k_);
    kinematic_model_.computeJacobian();
        
    fout.open("ik_solution.dat");
    if(fout) fout << "#iteration position_error_norm" << std::endl;
    //computeInitialGuess();
    computeSolution();
    
    solution = q_k_;
    int first_link = 0;
    int last_link = dofs_-1;
    kinematic_model_.setQ(solution);

    kinematic_model_.computeForwardKinematic(first_link, last_link);
    //kinematic_model_.computeJacobian();
   // kinematic_model_.getJacobian();
    std::cout << "trans in solveIk --> " << kinematic_model_.getTrans(first_link, last_link) << std::endl;
    fout.close();
}

void InverseKinematic::computeInitialGuess(){

    Eigen::MatrixXd jacobian;
    kinematic_model_.getJacobian(jacobian);
    Eigen::VectorXd q(dofs_);
    q = q_k_;

    int first_link = 0;
    int last_link = dofs_-1;
    kinematic_model_.computeForwardKinematic(first_link, last_link);
    Eigen::Vector3d pos_error = desired_pos_ - kinematic_model_.getTrans(first_link, last_link);

    double dpos_error_dt = 0;
    double pos_error_norm_prec = 0;

    bool is_ig_found = false;

    while(pos_error.norm()>EPSILON_ERROR && iteration < MAX_ITER_INITIAL_GUESS && !is_ig_found){

        Eigen::MatrixXd analitic_jacobian = jacobian.block(0,0,3,dofs_);
        Eigen::VectorXd dq = 0.8*analitic_jacobian.transpose()*pos_error;
        q += 0.1*dq;
        q_k_ = q_k_ + q;
        kinematic_model_.setQ(q);
        kinematic_model_.computeForwardKinematic(first_link, last_link);
        kinematic_model_.computeJacobian();
        kinematic_model_.getJacobian(jacobian);
        pos_error = desired_pos_ - kinematic_model_.getTrans(first_link, last_link);

        if(fout)
            fout << iteration << " " <<  pos_error.norm() << std::endl;
       
        iteration++;
        dpos_error_dt = (pos_error.norm()-pos_error_norm_prec)/0.001;
        pos_error_norm_prec = pos_error.norm();

        // if(std::abs(dpos_error_dt) < 0.0001){
        //     is_ig_found = true;
        // }
        
        // if(fout)
        //     fout << dpos_error_dt << std::endl;
          
    }

    std::cout << "q_k in IG --> " << q_k_ << std::endl;

}

void InverseKinematic::computeSolution(){
    for(int i=0; i<dofs_; i++){
        q_k_[i] = std::atan2(sin(q_k_[i]),cos(q_k_[i]));
        if(q_k_[i]<0) q_k_[i] += 2*M_PI;
    }
    Eigen::MatrixXd jacobian;
    kinematic_model_.setQ(q_k_);
    kinematic_model_.getJacobian(jacobian);

    int first_link = 0;
    int last_link = dofs_-1;
    kinematic_model_.computeForwardKinematic(first_link, last_link);
    Eigen::Vector3d pos_error = desired_pos_ - kinematic_model_.getTrans(first_link, last_link);

    while(pos_error.norm()>EPSILON_ERROR && iteration < MAX_ITER_INITIAL_SOL ){

        Eigen::MatrixXd analitic_jacobian = jacobian.block(0,0,3,dofs_);
        Eigen::MatrixXd pinv_jacobian = analitic_jacobian.completeOrthogonalDecomposition().pseudoInverse();
        Eigen::VectorXd dq = pinv_jacobian*pos_error;
        q_k_ +=  0.5*dq;
        kinematic_model_.setQ(q_k_);
        kinematic_model_.computeForwardKinematic(first_link, last_link);
        kinematic_model_.computeJacobian();
        kinematic_model_.getJacobian(jacobian);
        pos_error = desired_pos_ - kinematic_model_.getTrans(first_link, last_link);

        if(fout)
            fout << iteration << " " <<  pos_error.norm() << std::endl;
       
        iteration++;

          
    }

}

void InverseKinematic::initVariables(){

    dofs_ = robot_.getDofs();
    q_k_plus_one_.resize(dofs_);
    q_k_.resize(dofs_);

    desired_pos_.setZero();
    q_k_plus_one_.setZero();
    q_k_.setZero();
}

void InverseKinematic::setDesiredPos(Eigen::Vector3d &desired_pos){
    desired_pos_ = desired_pos;
}

void InverseKinematic::setQk(Eigen::VectorXd &q_k){
    q_k_ = q_k;
}

InverseKinematic::~InverseKinematic(){

}