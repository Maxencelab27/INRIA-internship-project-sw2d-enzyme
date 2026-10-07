#include <cmath>

#include "./swSource.h"
#include "../../../common/param.h"


double swSource::frictionChezy1D(const double& chezy_coeff, const double& hydraulic_radius, const double& u_norm, const double& dt)
/**
 * @brief Computes the factor by which the velocity (or unit discharge) is multiplied over one computational time step
 * @brief The formula is 1D, intended to be applied in the principal directions of the 2D friction tensor
 * @param[in] chezy_coeff: Chezy's friction coefficient
 * @param[in] hydraulic_radius: hydraulic radius (equal to h for swes, dip)
 * @param[in] u_norm: the norm of the flow velocity
 * @param[in] dt: current computational time step
 * @return The multiplication factor
 */
{
    double factor = 1.0;
    if ((hydraulic_radius > small_real) && (chezy_coeff > small_real))
    {
        double argument = - gravity * u_norm / chezy_coeff / chezy_coeff / hydraulic_radius;
        factor = std::exp(argument * dt);
    }
    else
    {
        factor = 0.0;
    }
    return factor;
}


double swSource::frictionManning1D(const double& manning_coeff, const double& hydraulic_radius, const double& u_norm, const double& dt)
/**
 * @brief Computes the factor by which the velocity (or unit discharge) is multiplied over one computational time step
 * @brief The formula is 1D, intended to be applied in the principal directions of the 2D friction tensor
 * @param[in] manning_coeff: Manning's friction coefficient
 * @param[in] hydraulic_radius: hydraulic radius (equal to h for swes, dip)
 * @param[in] u_norm: norm of the flow velocity
 * @param[in] dt: current computational time step
 * @return The multiplication factor
 */
{
    double factor = 1.0;
    if (hydraulic_radius > small_real)
    {
        double argument = - gravity * u_norm * manning_coeff * manning_coeff * std::exp(- 4.0 / 3.0 * std::log(hydraulic_radius));
        factor = std::exp(argument * dt);
    }
    else
    {
        factor = 0.0;
    }
    return factor;
}


double swSource::frictionStrickler1D(const double& strickler_coeff, const double& hydraulic_radius, const double& u_norm, const double& dt)
/**
 * @brief Computes the factor by which the velocity (or unit discharge) is multiplied over one computational time step
 * @brief The formula is 1D, intended to be applied in the principal directions of the 2D friction tensor
 * @param[in] strickler_coeff: Strickler's friction coefficient
 * @param[in] hydraulic_radius: hydraulic radius (equal to h for swes, dipp)
 * @param[in] u_norm: norm of the flow velocity
 * @param[in] dt: current computational time step
 * @return The multiplication factor
*/
{
    double factor = 1.0;
    if ((hydraulic_radius > small_real) && (strickler_coeff > small_real))
    {
        double argument = - gravity * u_norm / strickler_coeff / strickler_coeff * std::exp(-4.0 / 3.0 * std::log(hydraulic_radius));
        factor = std::exp(argument * dt);
    }
    else
    {
        factor = 0.0;
    }
    return factor;
}
