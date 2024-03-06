#ifndef ROBOT_H
#define ROBOT_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <forward_list>

#include <Eigen/Dense>
#include <json/json.h>

#include "link.h"

class Robot
{
    public:

        Robot();

        void buildRobotFromFile(std::string file_name);

        void getLinks(std::forward_list<Link> &links);

        int getDofs();

    private:

        int dofs_;
        // std::vector<Link> links_;
        std::forward_list<Link> links_;
      
};

#endif