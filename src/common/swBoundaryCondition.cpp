#include "swBoundaryCondition.h"

swBoundaryCondition::swBoundaryCondition(std::string type)
/**
 * @brief creates a boundary condition object (could be time series or rating curve)
 * @param[in] type of BC (q, h, z, ...)
 * @para Modified attribute:
 *  - bc_type
 */
{
    bc_type = type;
}

std::string swBoundaryCondition::getBCType()
/**
 * @brief returns the boundary condition type
 * @return bc_type (string) the boundary condition type
 */
{
    return bc_type;
}

double swBoundaryCondition::getBCVal()
/**
 * @brief returns the boundary condition value
 * @return bc_type (double) the boundary condition type
 */
{
    return bc_val;
}
