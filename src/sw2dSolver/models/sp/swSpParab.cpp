#include <algorithm>
#include <vector>

#include "./swSp.h"
#include "../../../common/param.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"
#include "../common/swParab.h"
#include "../common/swSource.h"


void swSp::parabComputeExplicit(void)
/**
 * @brief Sets up the matrix for the discretized diffusion equation
 * @return nothing
 * @par Modified attributes:
 *   - unit_discharge_x
 *   - unit_discharge_y
 */
{
    //---- Define edge diffusion coefficient as harmonic mean
    std::vector<double> diffusion_coefficient = {};
    for (size_t idx_edge = 0; idx_edge < edge_int_nb; idx_edge++)
    {
        size_t idx_celL = edge_int_celL_idx[idx_edge];
        size_t idx_celR = edge_int_celR_idx[idx_edge];
        double diff_L = momentum_diffusion_coefficient * water_depth[idx_celL] * phi[idx_celL];
        double diff_R = momentum_diffusion_coefficient * water_depth[idx_celR] * phi[idx_celR];
        double diff_edge = 0.0;
        if (std::max(diff_L, diff_R) > small_real)
        {
            diff_edge = 2.0 * diff_L * diff_R / (diff_L + diff_R);
        }
        diffusion_coefficient.push_back(diff_edge);
    }
    //==== Define edge diffusion coefficient as harmonic mean 

    //-- Diffusion of momentum_x
    parabComputeScalarExplicit(diffusion_coefficient, water_depth, unit_discharge_x);
    parabComputeScalarExplicit(diffusion_coefficient, water_depth, unit_discharge_y);
    //== Diffusion of momentum_x
}


/** @todo : adapt code for SP model (if needed) */
void swSp::parabComputeScalarExplicit(const std::vector<double>& diffusion_coefficient, const std::vector<double>& water_depth, std::vector<double>& scalar_cell)
/**
 * @brief Sets up the matrix for the discretized diffusion equation
 * @param[in] diffusion_coefficient: diffusion coefficient (defined edgewise, homogeneous in x-y)
 * @param[in] water_depth: water depth field (on cells)
 * @param[in,out] scalar_cell: any scalar field defined cellwise, to be modified by the diffusion process
 * @return nothing
 */
{

    //-- Define parabolic time step and number of substeps within one time step
    double dt_parab = t_max - t_0;
    for (size_t idx_edge = 0; idx_edge < edge_int_nb; idx_edge++)
    {
        size_t idx_celL = edge_int_celL_idx[idx_edge];
        size_t idx_celR = edge_int_celR_idx[idx_edge];
        dt_parab = std::min(dt_parab, 0.5 * cell_area[idx_celL] * phi[idx_celL] / diffusion_coefficient[idx_edge]);
        dt_parab = std::min(dt_parab, 0.5 * cell_area[idx_celR] * phi[idx_celR] / diffusion_coefficient[idx_edge]);
    }
    dt_parab = std::min(dt_parab, dt);
    size_t nb_parab_substeps = int(dt / dt_parab);
    dt_parab = dt / nb_parab_substeps;
    //== Define parabolic time step and number of substeps within one time step

    //-- Apply explicit scheme

    std::vector<double> flux(edge_int_nb, 0.0);
    //---- Loop over substeps
    for (size_t idx_substep = 0; idx_substep < nb_parab_substeps; idx_substep++)
    {
        //------ Loop over edges for fluxes
        for (size_t idx_edge = 0; idx_edge < edge_int_nb; idx_edge++)
        {
            size_t idx_celL = edge_int_celL_idx[idx_edge];
            size_t idx_celR = edge_int_celR_idx[idx_edge];
            flux[idx_edge] = diffusion_coefficient[idx_edge] * parab_edge_int_flux_factors[idx_edge] * 
                             (scalar_cell[idx_celL]/water_depth[idx_celL] - scalar_cell[idx_celR] / water_depth[idx_celR]);
        }
        //====== Loop over edges for fluxes
        //------ Loop over edges for balance
        for (size_t idx_edge = 0; idx_edge < edge_int_nb; idx_edge++)
        {
            size_t idx_celL = edge_int_celL_idx[idx_edge];
            size_t idx_celR = edge_int_celR_idx[idx_edge];
            scalar_cell[idx_celL] -= dt_parab * flux[idx_edge] / cell_area[idx_celL] / std::max(small_real, phi[idx_celL]);
            scalar_cell[idx_celR] += dt_parab * flux[idx_edge] / cell_area[idx_celR] / std::max(small_real, phi[idx_celR]);
        }
        //====== Loop over edges for balance
    }
    //---- Loop over substeps
    //== Apply explicit scheme
}