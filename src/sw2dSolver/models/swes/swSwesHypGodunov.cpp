#include <cmath>

#include "swSwes.h"
#include "../../../common/param.h"
#include "../../../common/swLog.h"

bool swSwes::hypGodunovComputeLRStates(void)
/**
 * @brief Computes the left and right states ( z, q, r, lambda+/- ) on every edge
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 *   - all interface-based Riemann states
 */
{
	//---- Boundary edges
    for (size_t idx_edge_bnd = 0 ; idx_edge_bnd < edge_bnd_nb ; idx_edge_bnd++)
    {
        double nx = edge_bnd_norm_x[idx_edge_bnd];
        double ny = edge_bnd_norm_y[idx_edge_bnd];
        size_t idx_cell = edge_bnd_cell_idx[idx_edge_bnd];
		double q = unit_discharge_x[idx_cell];
		double r = unit_discharge_y[idx_cell];
        h_riemann_bnd[idx_edge_bnd] = water_depth[idx_cell];
        z_riemann_bnd[idx_edge_bnd] = water_depth[idx_cell] + cell_z_b[idx_cell];
        q_riemann_bnd[idx_edge_bnd] = q * nx  + r * ny;
        r_riemann_bnd[idx_edge_bnd] = -q * ny + r * nx;
        double hR = h_riemann_bnd[idx_edge_bnd];
        double uR = q_riemann_bnd[idx_edge_bnd] / std::max(hR, small_real);
        double cR = std::sqrt(gravity * hR);
        lambda_plus_bnd[idx_edge_bnd] = lambdaPlus(uR, cR);
        lambda_minus_bnd[idx_edge_bnd] = lambdaMinus(uR, cR);
		lambda_max_bnd[idx_edge_bnd] = lambdaPlus(std::abs(uR) , cR);
	}
    //==== Boundary edges
	
    //---- Internal edges
    for (size_t idx_edge_int = 0 ; idx_edge_int < edge_int_nb ; idx_edge_int++)
    {
        double nx = edge_int_norm_x[idx_edge_int];
        double ny = edge_int_norm_y[idx_edge_int];
        
		// --- Left-hand side
        size_t id_celL = edge_int_celL_idx[idx_edge_int];
		double q = unit_discharge_x[id_celL];
		double r = unit_discharge_y[id_celL];
        h_riemann_int[idx_edge_int][0] = water_depth[id_celL];
        z_riemann_int[idx_edge_int][0] = water_depth[id_celL] + cell_z_b[id_celL];
        q_riemann_int[idx_edge_int][0] = q * nx  + r * ny;
        r_riemann_int[idx_edge_int][0] = -q * ny + r * nx;
        double hL = h_riemann_int[idx_edge_int][0];
        double uL = q_riemann_int[idx_edge_int][0] / std::max(hL, small_real);
        double cL = std::sqrt(gravity * hL);
        lambda_plus_int[idx_edge_int][0] = lambdaPlus(uL, cL);
        lambda_minus_int[idx_edge_int][0] = lambdaMinus(uL, cL);
		// === Left-hand side
		
		// --- Right-hand side
        size_t id_celR = edge_int_celR_idx[idx_edge_int];
		q = unit_discharge_x[id_celR];
		r = unit_discharge_y[id_celR];
        h_riemann_int[idx_edge_int][1] = water_depth[id_celR];
        z_riemann_int[idx_edge_int][1] = water_depth[id_celR] + cell_z_b[id_celR];
        q_riemann_int[idx_edge_int][1] = q * nx  + r * ny;
        r_riemann_int[idx_edge_int][1] = -q * ny + r * nx;
        double hR = h_riemann_int[idx_edge_int][1];
        double uR = q_riemann_int[idx_edge_int][1] / std::max(hR, small_real);
        double cR = std::sqrt(gravity * hR);
        lambda_plus_int[idx_edge_int][1] = lambdaPlus(uR, cR);
        lambda_minus_int[idx_edge_int][1] = lambdaMinus(uR, cR);
		// === Right-hand side
		lambda_max_int[idx_edge_int] = std::max(lambdaPlus(std::abs(uL) , cL) , lambdaPlus(std::abs(uR) , cR));
	}
    //==== internal edges
    return EXIT_SUCCESS;
}
