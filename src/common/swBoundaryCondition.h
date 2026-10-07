#pragma once

#include <string>

class swBoundaryCondition
/**
 * @brief The swBoundaryCondition class
 */
{
public:
    ///-- General Methods
    swBoundaryCondition(std::string type);
    ~swBoundaryCondition() = default;
    //== General Methods

    ///-- virtual methods
    virtual void computeValBC(const double& time) = 0;
    //== virtual methods

    ///-- getters
    double getBCVal();
    std::string getBCType();
    //== getters

protected:
    ///-- Attributes
    std::string bc_type;    /**< type of bc (h, q, ...) **/
    double bc_val;          /**< value of the parameter at the boundary (interpolated to the current time) **/
    //== Attributes
};
