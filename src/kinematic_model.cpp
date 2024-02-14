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

    R_ = Eigen::MatrixXd::Identity(3*dofs_,3*dofs_);
    trans_ = Eigen::VectorXd::Zero(3*dofs_);

    cos_theta = 0;
    sin_theta = 0;
    cos_alpha = 0;
    sin_alpha = 0;

    theta_prev_ = -1000;
    alpha_prev_ = -1000;

    //UR5 from https://www.universal-robots.com/articles/ur/application-installation/dh-parameters-for-calculations-of-kinematics-and-dynamics/
    //change also in computeForwardKinematic() if change robot
    // alpha_ << M_PI/2, 0, 0, M_PI/2, -M_PI/2, 0.0;
    // a_ << 0.0, -0.425, -0.39225, 0.0, 0.0, 0.0;
    // d_ << 0.089159, 0.0, 0.0, 0.10915, 0.09465, 0.0823;
}

void KinematicModel::computeForwardKinematic(const int start_link_idx, const int end_link_idx)
{
 
    Eigen::Matrix3d R; 
    Eigen::Vector3d trans;

    //double theta;
    // double cos_theta = cos(theta);
    // double sin_theta = sin(theta);
    

    for(int i=start_link_idx; i<end_link_idx; ++i){
        double theta = theta_(i) + q_(i);
        double a = a_(i);
        double d = d_(i);

        if(theta != theta_prev_)
        {
            sincosf(static_cast<float>(theta), &sin_theta, &cos_theta);
            theta_prev_ = theta;
        }
        if(alpha_(i) != alpha_prev_)
        {
            sincosf(static_cast<float>(alpha_(i)), &sin_alpha, &cos_alpha);
            alpha_prev_ = alpha_(i);
        }
        
        R << cos_theta, -sin_theta*cos_alpha, sin_theta*sin_alpha,
            sin_theta, cos_theta*cos_alpha , -cos_theta*sin_alpha,
            0           , sin_alpha              , cos_alpha;
        
        trans << a*cos_theta, a*sin_theta, d;

        trans_.segment(i*3,+3) += R_.block(i*3,i*3,3,3)*trans ;
        R_.block(i*3,i*3,3,3) *= R;
        
    }

}

void KinematicModel::computeJacobian(){

    Eigen::Vector3d pos_ee_absolute;
    pos_ee_absolute.setZero();
    computeForwardKinematic(0,dofs_-1);
    pos_ee_absolute = trans_;

    Eigen::Vector3d pos_ee_relative;
    pos_ee_relative.setZero();

    Eigen::Vector3d z_i_minus_one;
    z_i_minus_one.setZero();

    Eigen::Vector3d z_i_minus_one_const;
    z_i_minus_one_const << 0,0,1;

    for(int i=0; i<dofs_; i++){
        computeForwardKinematic(0,i); 
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

Eigen::Vector3d KinematicModel::getTrans(const int start_link_idx, const int end_link_idx){
    return trans_.segment(start_link_idx*3,end_link_idx*3-start_link_idx*3+3);
}

Eigen::Matrix3d KinematicModel::getR(const int start_link_idx, const int end_link_idx){
    return R_.block(start_link_idx*3,start_link_idx*3,end_link_idx*3-start_link_idx*3+3,end_link_idx*3-start_link_idx*3+3);
}

void KinematicModel::getJacobian(Eigen::MatrixXd &jacobian){
    jacobian = jacobian_;
}

KinematicModel::~KinematicModel(){

}