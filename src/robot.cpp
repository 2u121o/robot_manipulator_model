#include "robot.h"


Robot::Robot(){

    std::cout << "Robot constructed! Fill the structure." << std::endl;

}

void Robot::buildRobotFromFile(std::string file_name){

    std::ifstream parameters_file(file_name);

    Json::Value root;
    Json::Reader reader;
    bool is_parsing_succesful = reader.parse(parameters_file, root);

    if(!is_parsing_succesful){
        //TODO: lauch an exception
        std::cout << "Failed to parse JSON: " << reader.getFormattedErrorMessages() << std::endl;
    }

    std::cout << "File parsed successfully" << std::endl;

    std::string robot_name = root["robot"]["name"].asString();
    dofs_ = root["robot"]["dofs"].asInt();

    links_.resize(dofs_);

    Json::Value links = root["robot"]["links"];
    for(int link_number=0; link_number<dofs_; link_number++){

        DHParams dh_params;
        Json::Value dh_value = links[link_number]["dh"];
        dh_params.theta = dh_value["theta"].asDouble();
        dh_params.d = dh_value["d"].asDouble();
        dh_params.a = dh_value["a"].asDouble();
        dh_params.alpha = dh_value["alpha"].asDouble();

        DynamicParameters dynamic_parameters;
        Json::Value dynamic_parameters_value = links[link_number]["dynamic_params"];
        dynamic_parameters.mass = dynamic_parameters_value["mass"].asDouble();

        Eigen::Vector3d com;
        com(0) = dynamic_parameters_value["com"][0].asDouble();
        com(1) = dynamic_parameters_value["com"][1].asDouble();
        com(2) = dynamic_parameters_value["com"][2].asDouble();
        dynamic_parameters.com = com;

        Eigen::VectorXd inertia;
        const int NUM_INERTIA_PARAMS = 6;
        inertia.resize(NUM_INERTIA_PARAMS);
        Json::Value inertia_value = dynamic_parameters_value["inertia"];
        inertia(0) = inertia_value["Ixx"].asDouble();
        inertia(1) = inertia_value["Ixy"].asDouble();
        inertia(2) = inertia_value["Ixz"].asDouble();
        inertia(3) = inertia_value["Iyy"].asDouble();
        inertia(4) = inertia_value["Iyz"].asDouble();
        inertia(5) = inertia_value["Izz"].asDouble();
        dynamic_parameters.inertia = inertia;

        Link link(link_number, dh_params, dynamic_parameters);
        links_.at(link_number) = link;
    } 

    parameters_file.close();

}

void Robot::getLinks(std::vector<Link> &links){
    links = links_;
}

int Robot::getDofs(){
    return dofs_;
}

Robot::~Robot(){

}

        
