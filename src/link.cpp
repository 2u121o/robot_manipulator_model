#include "link.h"

Link::Link()
{

}

Link::Link(const int link_number, const DHParams &dh_params, const DynamicParameters &dynamic_parameters):
link_number_(link_number),
dh_params_(dh_params),
dynamic_parameters_(dynamic_parameters)
{

}

void Link::setLinkNumber(const int link_number)
{
    link_number_ = link_number;
}

void Link::setDHParams(const DHParams &dh_params)
{
    dh_params_ = dh_params;
}

void Link::setDynamicParameters(const DynamicParameters &dynamic_parameters)
{
    dynamic_parameters_ = dynamic_parameters;
}

int Link::getLinkNumber()
{
    return link_number_;
}

void Link::getDHParams(DHParams &dh_params)
{
    dh_params = dh_params_;
}

void Link::getDynamicParameters(DynamicParameters &dynamic_parameters)
{
    dynamic_parameters = dynamic_parameters_;
}
