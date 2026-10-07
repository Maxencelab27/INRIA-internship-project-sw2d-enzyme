#include <algorithm>
#include <cmath>

#include "./swSource.h"
#include "../../../common/param.h"
#include "../../../common/swLog.h"


void swSource::windConvertFromPolarToCartesian(const std::vector<double>& norm, const std::vector<double>& angle_degrees,
                                               std::vector<double>& component_x, std::vector<double>& component_y)
/**
 * @brief Converts a table of wind time series in the form (norm, angle_from_North) to tables of the x- and y-components
 * @param[in] norm[n]: norm of the wind velocity for the n-th time record
 * @param[in] angle_degrees[n]: angle (degrees) from North of the wind velocity for the n-th time record
 * @param[out] component_x[n]: x-velocity at the n-th time record
 * @param[out] component_y[n]: y-velocity at the n-th time record
 */
{
    //-- Ensure consistency
    if (norm.size() != angle_degrees.size())
    {
        logError() << "Wind norm and angle vectors do not have the same size - cannot proceed with conversion";
        logFatal(__FILE__, __func__, __LINE__);
    }
    size_t n_records = norm.size();
    component_x.resize(n_records);
    component_y.resize(n_records);
    //== Ensure consistency

    //-- Perform conversion (NB. Angles in degrees, counted clockwise from North)
    for (size_t n = 0; n < n_records; n++)
    {
        double angle_radians = atan(1.0) / 45. * (90. - angle_degrees[n]);
        component_x[n] = -norm[n] * cos(angle_radians);
        component_y[n] = -norm[n] * sin(angle_radians);
    }
    //== Perform conversion (NB. Angles in degrees, counted clockwise from North)
}
