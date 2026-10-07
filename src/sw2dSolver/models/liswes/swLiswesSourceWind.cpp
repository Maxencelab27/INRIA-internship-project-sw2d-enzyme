#include <cmath>
#include <vector>

#include "./swLiswes.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"
#include "../common/swSource.h"
#include "../../../common/param.h"


void swLiswes::sourceWind()
/**
 * @brief advances the solution in time for the wind stress step
 * @return nothing
 * @par Modified attributes:
 *   - unit_discharge_x
 *   - unit_discharge_y
 */
{
    //-- Interpolate wind time series
    size_t wind_station_nb = wind_velocity_x_table.size();
    std::vector<double> stationwise_wind_velocity_x_interpolated = {};
    stationwise_wind_velocity_x_interpolated.resize(wind_station_nb);
    std::vector<double> stationwise_wind_velocity_y_interpolated = {};
    stationwise_wind_velocity_y_interpolated.resize(wind_station_nb);

    for (size_t idx_station = 0; idx_station < wind_station_nb; idx_station++)
    {
        // VGt 22.06.2025 - Using t + dt/2 yields 2nd-order time accuracy
        stationwise_wind_velocity_x_interpolated[idx_station] = interpolateTimeSeries(wind_times_vector, wind_velocity_x_table[idx_station], t +0.5 * dt);
        stationwise_wind_velocity_y_interpolated[idx_station] = interpolateTimeSeries(wind_times_vector, wind_velocity_y_table[idx_station], t + 0.5 * dt);
    }
    //== Interpolate wind time series

    if (model_wind == "smithbanke")
    {
        sourceWindSmithBanke(stationwise_wind_velocity_x_interpolated, stationwise_wind_velocity_y_interpolated);
    }
    else
    {
        logError() << "Wind model unknown : " << model_wind;
        logFatal(__FILE__, __func__, __LINE__);
    }
}


void swLiswes::sourceWindSmithBanke(const std::vector<double>& stationwise_wind_velocity_x_interpolated, 
                                  const std::vector<double>& stationwise_wind_velocity_y_interpolated)
/**
 * @brief Advances the solution in time for wind stress step (Smith & Banke's model)
 * @brief Modifies the unit_discharge_x and unit_discharge_y vectors
 * @param[in] stationwise_wind_velocity_x_interpolated: x-velocity that has been interpolated from the tables at the current time
 * @param[in] stationwise_wind_velocity_y_interpolated: y-velocity that has been interpolated from the tables at the current time
 * @return nothing
 * @par Modified attributes:
 *   - unit_discharge_x
 *   - unit_discharge_y
 */
{

    //---- Loop on cells to compute the source term
    for (size_t idx_cell = 0 ; idx_cell < cell_nb ; idx_cell++)
    {
        size_t idx_station = wind_station_code[idx_cell] - 1; // NB: code start at 1, table starts with index 0 !
        double x_wind_velocity = stationwise_wind_velocity_x_interpolated[idx_station];
        double y_wind_velocity = stationwise_wind_velocity_y_interpolated[idx_station];
        double h0 = water_depth[idx_cell];
        double x_velocity = unit_discharge_x[idx_cell] / h0;
        double y_velocity = unit_discharge_y[idx_cell] / h0;
        double x_relative_velocity = x_wind_velocity - x_velocity;
        double y_relative_velocity = y_wind_velocity - y_velocity;
        double wind_norm = std::sqrt(x_velocity * x_velocity + y_velocity * y_velocity);
        double wind_relative_norm = std::sqrt(x_relative_velocity * x_relative_velocity + y_relative_velocity * y_relative_velocity);

        if (h0 > small_real) {
            double CsD = 1e-3 * (0.63 + 0.066 * wind_norm);
            double tau_x = rho_air * CsD * wind_relative_norm * x_relative_velocity;
            double mcoef = std::min(h0, h_min) / std::max(h_min, small_real);
            double q_new_x = unit_discharge_x[idx_cell] + tau_x / rho_water * dt * mcoef;
            double tau_y = rho_air * CsD * wind_relative_norm * y_relative_velocity;
            double q_new_y = unit_discharge_y[idx_cell] + tau_y / rho_water * dt * mcoef;
            unit_discharge_x[idx_cell] = q_new_x;
            unit_discharge_y[idx_cell] = q_new_y;
        }
    }
}