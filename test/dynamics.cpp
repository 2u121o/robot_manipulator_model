#include <iostream>

#include <urdf/model.h>
#include <kdl/chain.hpp>
#include <kdl/chaindynparam.hpp>
#include <kdl/frames.hpp>
#include <kdl/tree.hpp>
#include <kdl_parser/kdl_parser.hpp>

int main()
{

    urdf::Model urdf_model;
    urdf_model.initFile("/home/dario/Personal/Workspace/robot_manipulator_model/urdfs/UR5_urdf.urdf");     
    KDL::Tree robot_tree;
    kdl_parser::treeFromUrdfModel(urdf_model, robot_tree);

    KDL::Chain robot_chain;
    robot_tree.getChain("base", "link5", robot_chain);

    KDL::Vector gravity(0.0, 0.0, -9.81);

    KDL::JntArray q(robot_chain.getNrOfJoints());
    q(0) = 1.0;
    KDL::ChainDynParam dyn_param(robot_chain, gravity);
    KDL::JntSpaceInertiaMatrix mass_matrix(robot_chain.getNrOfJoints());
    dyn_param.JntToMass(q, mass_matrix);

    KDL::JntArray q_dot(robot_chain.getNrOfJoints());

    KDL::JntArray coriolis_centrifugal_forces(robot_chain.getNrOfJoints());
    dyn_param.JntToCoriolis(q, q_dot, coriolis_centrifugal_forces);

    KDL::JntArray gravity_forces(robot_chain.getNrOfJoints());
    dyn_param.JntToGravity(q, gravity_forces);
    std::cout << "ci sono gravity_forces.rows() " << gravity_forces.rows() <<std::endl;
    for(int i=0; i<gravity_forces.rows(); ++i)
    {
        std::cout << "Joint " << i << ": " << gravity_forces(i) << " N" << std::endl;
    }


}