#ifndef ROBOT_H
#define ROBOT_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <Eigen/Dense>
#include <json/json.h>

#include "link.h"

class Robot
{
    public:

        Robot();
        ~Robot();

        void buildRobotFromFile(std::string file_name);

        void getLinks(std::vector<Link> &links);

        int getDofs();

    private:

        int dofs_;
        std::vector<Link> links_;
      
};

#endif