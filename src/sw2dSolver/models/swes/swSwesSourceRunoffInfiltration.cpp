#include <cmath>
#include <vector>

#include "./swSwes.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"
#include "../common/swSource.h"
#include "../../../common/param.h"


void swSwes::sourceInfiltration(void)
/**
 * @brief advances the solution in time for the infiltration step
 * @return nothing
 * @par Modified attributes:
 *   - water_depth
 *   - unit_discharge_x
 *   - unit_discharge_y
 */
{
    if (model_infiltration == "h_proportional")
    {
        logError() << "swSwes::sourceInfiltration - model_infiltration 'h_proportionnal' available only for industrial use (CEREG) ";
        logError() << "Please set another value for the model_infiltration parameter. ";
        logFatal(__FILE__, __func__, __LINE__);
    }
    else if (model_infiltration == "horton")
    {
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            double h_0 = water_depth[idx_cell];
            double dh = - infiltration_param[idx_cell]  * dt;
            water_depth[idx_cell] = std::max(0.0, h_0 + dh);
            double ratio = water_depth[idx_cell] / std::max(h_0, small_real);
            unit_discharge_x[idx_cell] *= ratio;
            unit_discharge_y[idx_cell] *= ratio;
        }
    }
}


void swSwes::sourceRunoff(void)
/**
 * @brief advances the solution in time for the runoff generation step
 * @return nothing
 * @par Modified attributes:
 *   - water_depth
 *   - unit_discharge_x
 *   - unit_discharge_y
 */
{
    //-- Interpolate precipitation time series
    size_t precipitation_station_nb = stationwise_precipitation_rate_table.size();
    std::vector<double> stationwise_precipitation_rate_interpolated = {};
    stationwise_precipitation_rate_interpolated.resize(precipitation_station_nb);

    for (size_t idx_station = 0; idx_station < precipitation_station_nb; idx_station++)
    {
        // VGt 22.06.2025 - Using t + dt/2 yields 2nd-order time accuracy
        stationwise_precipitation_rate_interpolated[idx_station] = interpolateTimeSeries(precipitation_times_vector,
                                                                                         stationwise_precipitation_rate_table[idx_station],
                                                                                         t + 0.5 * dt);
    }
    //== Interpolate precipitation time series

    //-- Apply Runoff (including negative P rate)
    if (model_runoff == "full_rainfall")
    {
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            double h0 = water_depth[idx_cell]; // Initial water depth
            double dp = dt * stationwise_precipitation_rate_interpolated[precipitation_station_code[idx_cell] - 1]; // net rainfall depth over time step
            if (dp > 0.) // if rainfall, only 1 - inf. coeff contributes to runoff
            {
                water_depth[idx_cell] += dp;
            }
            else // if evaporation, then 100% of water on surface returns to the atmosphere, momentum must be updated
            {
                water_depth[idx_cell] = std::max(0.0, h0 + dp);
                double ratio = water_depth[idx_cell] / std::max(h0, small_real);
                unit_discharge_x[idx_cell] *= ratio;
                unit_discharge_y[idx_cell] *= ratio;
            }
        }
    }
    else if (model_runoff == "infiltration_coefficient")
    {
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            double h0 = water_depth[idx_cell]; // Initial water depth
            double dp = dt * stationwise_precipitation_rate_interpolated[precipitation_station_code[idx_cell] -1]; // net rainfall depth over time step
            double dh = 0.0;
            if (dp > 0.) // if rainfall, only 1 - inf. coeff contributes to runoff
            {
                dh = (1.0 - runoff_param[idx_cell]) * dp;
                water_depth[idx_cell] += dh;
            }
            else // if evaporation, then 100% of water on surface returns to the atmosphere, momentum must be updated
            {
                dh = dp;
                water_depth[idx_cell] = std::max(0.0, h0 + dh);
                double ratio = water_depth[idx_cell] / std::max(h0, small_real);
                unit_discharge_x[idx_cell] *= ratio;
                unit_discharge_y[idx_cell] *= ratio;
            }
        }
    }
    else if (model_runoff == "runoff_coefficient")
    {
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            double h0 = water_depth[idx_cell]; // Initial water depth
            double dp = dt * stationwise_precipitation_rate_interpolated[precipitation_station_code[idx_cell] - 1]; // net rainfall depth over time step
            double dh = 0.0;
            if (dp > 0.) // if rainfall, only runoff_coeff contributes to runoff
            {
                dh = runoff_param[idx_cell] * dp;
                water_depth[idx_cell] += dh;
            }
            else // if evaporation, then 100% of water on surface returns to the atmosphere, momentum must be updated
            {
                dh = dp;
                water_depth[idx_cell] = std::max(0.0, h0 + dh);
                double ratio = water_depth[idx_cell] / std::max(h0, small_real);
                unit_discharge_x[idx_cell] *= ratio;
                unit_discharge_y[idx_cell] *= ratio;
            }
        }
    }
    else if (model_runoff == "scs")
    {
        logError() << "swSwes::sourceRunoff - model_runoff 'scs' available only for industrial use (CEREG) ";
        logError() << "Please set another value for the model_runoff parameter. ";
        logFatal(__FILE__, __func__, __LINE__);
    }
    else
    {
        logError() << "swSwes::sourceRunoff - model_runoff not handled = " << model_runoff;
        logFatal(__FILE__, __func__, __LINE__);
    }
    //== Apply Runoff (including negative P rate)
}