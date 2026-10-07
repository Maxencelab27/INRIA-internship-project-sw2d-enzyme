#include "swBoundaryConditionTimeSeries.h"

#include "../../../common/swLog.h"


swBoundaryConditionTimeSeries::swBoundaryConditionTimeSeries(std::string type) : swBoundaryCondition(type)
/**
 * @brief constructor of the swBoundaryConditionTimeSeries class
 * @param[in] type: type of the bc: q, h, z, fr, wall
 * @param[in] param: pointer to the swParamSolver object
 */
{
    // do nothing on purpose
}


void swBoundaryConditionTimeSeries::initialise(const std::vector<double> &time, const std::vector<double> &val)
/**
 * @brief initialise the boundary and define the values for the boundary flux computation
 * @param[in] time: vector containing the time values for the current bc
 * @param[in] val: vector containing the bc values for the current bc
 */
{
    if (val.size() != time.size())
    {
        logError() << "Error while building bc";
        logError() << "  " << time.size() << "time values";
        logError() << "  " << val.size() << "variable values";
        logError() << "  Should be the same";
        logFatal(__FILE__, __func__, __LINE__);
    }

    if (val.size() < 1)
    {
        logError() << "Error while building bc";
        logError() << "  At least ont value have to be defined for the boundary condition";
        logError() << "  " << val.size() << "values read";
        logFatal(__FILE__, __func__, __LINE__);
    }

    time_values = time;
    bc_values = val;

    v_inf = bc_values[0];
    t_inf = time_values[0];
    if (time_values.size() > 1)
    {
        v_sup = bc_values[1];
        t_sup = time_values[1];
    }
    else
    {
        v_sup = v_inf;
        t_sup = t_inf;
        at_end = true;
    }
    rank = 1;
}


void swBoundaryConditionTimeSeries::computeValBC(const double& time)
/**
 * @brief interpolate the value at the bc for the current time
 * @param[in] time: current simulation time
 * @par Modified attributes:
 *  - bc_val 
 */
{
    while ((time > t_sup) && (at_end == false))
    {
        v_inf = v_sup;
        t_inf = t_sup;
        t_sup = time_values[rank];
        v_sup = bc_values[rank];
        ++rank;

        if (rank == time_values.size())
        {
            at_end = true;
        }
    }

    if (time > t_sup)
    {
        bc_val = v_sup;
    }
    else if (time < t_inf)
    {
        bc_val = v_inf;
    }
    else
    {
        double delta_t = (time - t_inf) / (t_sup - t_inf);
        bc_val = v_inf + (v_sup - v_inf) * delta_t;
    }
}
