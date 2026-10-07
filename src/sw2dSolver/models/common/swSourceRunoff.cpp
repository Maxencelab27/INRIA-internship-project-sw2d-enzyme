#include <algorithm>
#include <cmath>

#include "./swSource.h"
#include "../../../common/param.h"
#include "../../../common/swLog.h"


double swSource::runoffScs(const double& dt, const double& i_a, const double& s, const double& p, double& h_a, double& h_s)
/**
 * @brief Computes runoff generation using the 2-tank, conceptual model version of the SCS model
 * @param[in] dt: computational time step
 * @param[in] i_a: capacity of the initial accumulation compartment
 * @param[in] s: capacity of the retention compartment
 * @param[in] p: precipitation rate
 * @param[out] h_a: water depth in the i_a compartment
 * @param[out] h_s: water depth in the s compartment
 * @return r: runoff intensity
 */
{
    //-- Store initial values
    double h_a_0 = h_a;
    double h_s_0 = h_s;
    //== Store initial values

    //-- Compute next h_a value and flux q leaving i_a
    h_a = std::max(0.0, h_a_0 + dt * p);
    double q = std::max(0.0, h_a - i_a) / dt;
    h_a = std::min(h_a, i_a);
    //== Compute next h_a value and flux q leaving i_a

    //-- Compute next h_s value and flux r leaving system
    double x_0 = 1.0 - h_s / std::max(s, small_real);
    double x = 1.0 / (1.0 / std::max(small_real, x_0) + q * dt / std::max(s, small_real) );
    h_s = (1.0 - x) * s;
    double r = q + (h_s_0 - h_s) / dt;
    //== Compute next h_s value and flux r leaving system
    return r;
}
