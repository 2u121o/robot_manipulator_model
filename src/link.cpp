#include "link.h"

Link::Link(std::string name, double mass,
     Eigen::VectorXd inertia, Eigen::Vector3d cog, Eigen::Vector3d dh)
    :name_(name), inertia_(inertia), cog_(cog), dh_(dh)
{
    if(mass>0){
        mass_ = mass;
        std::cout << "Link initialized" << std::endl;
    }
    else
        std::cerr << "Mass must be positive" << std::endl;


}

Link::Link(std::string name)
    :name_(name)
{
    std::cout << "Link initialized just with name" << std::endl;
}

void Link::setName(std::string name){
    name_ = name;
}
void Link::setMass(double mass){
    if(mass > 0)
        mass_ = mass;
    else
        std::cout << "Mass must be positive" << std::endl;
}
void Link::setInertia(Eigen::VectorXd inertia){
        inertia_ = inertia;
}
void Link::setCog(Eigen::Vector3d cog){
    cog_ = cog;
}
void Link::setDh(Eigen::Vector3d dh){
    dh_ = dh;
}

std::string Link::getName(){
    return name_;
}
double Link::getMass(){
    return mass_;
}
Eigen::VectorXd Link::getInertia(){
    return inertia_;
}
Eigen::Vector3d Link::getCog(){
    return cog_;
}
Eigen::Vector3d Link::getDh(){
    return dh_;
}
