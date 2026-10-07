#include <algorithm>
#include <vector>

#include "./swDdp.h"
#include "../../../common/param.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"
#include "../common/swParab.h"
#include "../common/swSource.h"


/** @todo: adapt code */
void swDdp::parabComputeExplicit(void)
/**
 * @brief Sets up the matrix for the discretized diffusion equation
 * @return nothing
 * @par Modified attributes:
 *   - unit_discharge_x
 *   - unit_discharge_y
 */
{
}


/** @todo: adapt code */
void swDdp::parabComputeScalarExplicit(const std::vector<double>& diffusion_coefficient, const std::vector<double>& water_depth, std::vector<double>& scalar_cell)
/**
 * @brief Sets up the matrix for the discretized diffusion equation
 * @param[in] diffusion_coefficient: diffusion coefficient (defined edgewise, homogeneous in x-y)
 * @param[in] water_depth: water depth field (on cells)
 * @param[in,out] scalar_cell: any scalar field defined cellwise, to be modified by the diffusion process
 * @return nothing
 */
{
}