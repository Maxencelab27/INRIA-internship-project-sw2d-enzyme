#pragma once

#include <string>
#include <vector>
#include "../../../common/swBoundaryCondition.h"

class swBoundaryConditionTimeSeries : public swBoundaryCondition
/**
 * @brief The swBoundaryConditionTimeSeries class
 */
{
public:
    //-- General Methods
    swBoundaryConditionTimeSeries(std::string type);

    void initialise(const std::vector<double> &time, const std::vector<double> &val);
    void computeValBC(const double& time) final;
    //== General Methods

    std::vector<double> time_values;    /**< list of the time values for interpolation **/
    std::vector<double> bc_values;      /**< list of the values for interpolation **/

protected:
    //-- Attributes
    bool at_end = false;    /**< at the end of the list (no interpolation) **/
    size_t rank = 0;        /**< current rank in the rating curve**/
    double t_inf = 0.;      /**< previous value of t (for interpolation) **/
    double t_sup = 0.;      /**< next value of t (for interpolation) **/
    double v_inf = 0.;      /**< previous value (for interpolation) **/
    double v_sup = 0.;      /**< next value (for interpolation) **/
    //== Attributes 
};
