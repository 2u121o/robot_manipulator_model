#include "kinematic_model.h"

KinematicModel::KinematicModel(){

}

KinematicModel::KinematicModel(Robot robot):robot_(robot){

    initVariables();
    std::cout << "Kinematic model constructed! " << std::endl;
}

void KinematicModel::initVariables(){

    std::vector<Link> links;
    robot_.getLinks(links);

    dofs_=  robot_.getDofs();

    q_.resize(dofs_);
    jacobian_.resize(dofs_, dofs_);
    theta_.resize(dofs_);
    d_.resize(dofs_);
    a_.resize(dofs_);
    alpha_.resize(dofs_);

    q_.setZero();
    cartesina_pos_.setZero();
    jacobian_.setZero();

    for(int i=0; i<dofs_; i++){

        Link link = links.at(i);

        DHParams dh_params;
        link.getDHParams(dh_params);

        theta_(i) = dh_params.theta;
        d_(i) = dh_params.d;
        a_(i) = dh_params.a;
        alpha_(i) = dh_params.alpha;

    }

    //UR5 from https://www.universal-robots.com/articles/ur/application-installation/dh-parameters-for-calculations-of-kinematics-and-dynamics/
    //change also in computeForwardKinematic() if change robot
    // alpha_ << M_PI/2, 0, 0, M_PI/2, -M_PI/2, 0.0;
    // a_ << 0.0, -0.425, -0.39225, 0.0, 0.0, 0.0;
    // d_ << 0.089159, 0.0, 0.0, 0.10915, 0.09465, 0.0823;
}

void KinematicModel::computeForwardKinematic(const std::vector<int> idx_links){
 
    R_.setIdentity();
    trans_.setZero();

    Eigen::Matrix3d R; 
    R.setZero();
    Eigen::Vector3d trans;
    trans.setZero();

    //double theta;

    for(int i=idx_links.at(0); i<idx_links.at(1); i++){

        double theta = theta_(i);
        double alpha = alpha_(i);
        double a = a_(i);
        double d = d_(i);

        theta += q_[i];

        R << cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha),
            sin(theta), cos(theta)*cos(alpha) , -cos(theta)*sin(alpha),
            0           , sin(alpha)              , cos(alpha);

        trans << a*cos(theta), a*sin(theta), d;

        trans_ = R_*trans + trans_;
        R_ = R_*R;
    }

}

void KinematicModel::computeJacobian(){

    Eigen::Vector3d pos_ee_absolute;
    pos_ee_absolute.setZero();
    std::vector<int> link_origins = {0,dofs_-1};
    computeForwardKinematic(link_origins);
    pos_ee_absolute = trans_;

    Eigen::Vector3d pos_ee_relative;
    pos_ee_relative.setZero();

    Eigen::Vector3d z_i_minus_one;
    z_i_minus_one.setZero();

    Eigen::Vector3d z_i_minus_one_const;
    z_i_minus_one_const << 0,0,1;

    for(int i=0; i<dofs_; i++){
        std::vector<int> link_origins = {0,i};
        computeForwardKinematic(link_origins); //lerrore e probabilmente qui verifica questa parte 
        z_i_minus_one = R_*z_i_minus_one_const;
        pos_ee_relative = pos_ee_absolute - trans_;
        jacobian_.block(0,i,3,1) = z_i_minus_one.cross(pos_ee_relative);
        jacobian_.block(3,i,3,1) = z_i_minus_one;
    }
}

void KinematicModel::setQ(const Eigen::VectorXd &q){
    q_ = q;
}

Eigen::VectorXd KinematicModel::getQ(){
    return q_;
}

Eigen::Vector3d KinematicModel::getTrans(){
    return trans_;
}

Eigen::Matrix3d KinematicModel::getR(){
    return R_;
}

void KinematicModel::getJacobian(Eigen::MatrixXd &jacobian){
    jacobian = jacobian_;
}

KinematicModel::~KinematicModel(){

}