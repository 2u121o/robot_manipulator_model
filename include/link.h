#ifndef LINK_H
#define LINK_H

#include <iostream>

#include <Eigen/Dense>

#include "utils.h"

template <typename T>
class Link
{

    public:

        Link();
        Link(const int link_number, const DHParams<T> &dh_params, const DynamicParameters<T> &dynamic_parameters);

        void setLinkNumber(const int link_number);
        void setDHParams(const DHParams<T> &dh_params);
        void setDynamicParameters(const DynamicParameters<T> &dynamic_parameters);

        int getLinkNumber();
        void getDHParams(DHParams<T> &dh_params);
        void getDynamicParameters(DynamicParameters<T> &dynamic_parameters);
 
    private:

        int link_number_;
        DHParams<T> dh_params_;
        DynamicParameters<T> dynamic_parameters_;
    
};

template <typename T>
Link<T>::Link()
{

}

template <typename T>
Link<T>::Link(const int link_number, const DHParams<T> &dh_params, const DynamicParameters<T> &dynamic_parameters):
link_number_(link_number),
dh_params_(dh_params),
dynamic_parameters_(dynamic_parameters)
{

}

template <typename T>
void Link<T>::setLinkNumber(const int link_number)
{
    link_number_ = link_number;
}

template <typename T>
void Link<T>::setDHParams(const DHParams<T> &dh_params)
{
    dh_params_ = dh_params;
}

template <typename T>
void Link<T>::setDynamicParameters(const DynamicParameters<T> &dynamic_parameters)
{
    dynamic_parameters_ = dynamic_parameters;
}

template <typename T>
int Link<T>::getLinkNumber()
{
    return link_number_;
}

template <typename T>
void Link<T>::getDHParams(DHParams<T> &dh_params)
{
    dh_params = dh_params_;
}

template <typename T>
void Link<T>::getDynamicParameters(DynamicParameters<T> &dynamic_parameters)
{
    dynamic_parameters = dynamic_parameters_;
}

#endif