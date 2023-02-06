#ifndef LINK_H
#define LINK_H

#include <iostream>

#include <Eigen/Dense>

#include "utils.h"

class Link{

    public:

        Link();
        Link(int link_number, DHParams dh_params, DynamicParameters dynamic_parameters);
        ~Link();

        void setLinkNumber(const int link_number);
        void setDHParams(const DHParams &dh_params);
        void setDynamicParameters(const DynamicParameters dynamic_parameters);

        int getLinkNumber();
        void getDHParams(DHParams &dh_params);
        void getDynamicParameters(DynamicParameters &dynamic_parameters);
 
    private:

        int link_number_;
        DHParams dh_params_;
        DynamicParameters dynamic_parameters_;
        

};

#endif