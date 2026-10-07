#include <cmath>

#include "./swLiswes.h"
#include "../../../common/param.h"
#include "../../../common/swLog.h"

bool swLiswes::hypGodunovComputeLRStates(void)
/**
 * @brief Computes the left and right states ( z, q, r, lambda+/- ) on every edge
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 *   - all interface-based Riemann states
 */
{
	//---- Boundary edges
    for (size_t edge_idx = 0 ; edge_idx < edge_bnd_nb ; edge_idx++)
    {
        double nx = edge_bnd_norm_x[edge_idx];
        double ny = edge_bnd_norm_y[edge_idx];
        size_t id_cell = edge_bnd_cell_idx[edge_idx];
		double q = unit_discharge_x[id_cell];
		double r = unit_discharge_y[id_cell];
        h_riemann_bnd[edge_idx] = water_depth[id_cell];
        z_riemann_bnd[edge_idx] = water_depth[id_cell] + cell_z_b[id_cell];
        q_riemann_bnd[edge_idx] = q * nx  + r * ny;
        r_riemann_bnd[edge_idx] = -q * ny + r * nx;
        double hR = h_riemann_bnd[edge_idx];
        double cR = std::sqrt(gravity * hR);
        lambda_plus_bnd[edge_idx] = lambdaPlus(cR);
        lambda_minus_bnd[edge_idx] = lambdaMinus(cR);
		lambda_max_bnd[edge_idx] = lambdaPlus(cR);
	}
    //==== Boundary edges
	
    //---- Internal edges
    for (size_t edge_idx = 0 ; edge_idx < edge_int_nb ; edge_idx++)
    {
        double nx = edge_int_norm_x[edge_idx];
        double ny = edge_int_norm_y[edge_idx];
        
		// --- Left-hand side
        size_t id_celL = edge_int_celL_idx[edge_idx];
		double q = unit_discharge_x[id_celL];
		double r = unit_discharge_y[id_celL];
        h_riemann_int[edge_idx][0] = water_depth[id_celL];
        z_riemann_int[edge_idx][0] = water_depth[id_celL] + cell_z_b[id_celL];
        q_riemann_int[edge_idx][0] = q * nx  + r * ny;
        r_riemann_int[edge_idx][0] = -q * ny + r * nx;
        double hL = h_riemann_int[edge_idx][0];
        double cL = std::sqrt(gravity * hL);
        lambda_plus_int[edge_idx][0] = lambdaPlus(cL);
        lambda_minus_int[edge_idx][0] = lambdaMinus(cL);
		// === Left-hand side
		
		// --- Right-hand side
        size_t id_celR = edge_int_celR_idx[edge_idx];
		q = unit_discharge_x[id_celR];
		r = unit_discharge_y[id_celR];
        h_riemann_int[edge_idx][1] = water_depth[id_celR];
        z_riemann_int[edge_idx][1] = water_depth[id_celR] + cell_z_b[id_celR];
        q_riemann_int[edge_idx][1] = q * nx  + r * ny;
        r_riemann_int[edge_idx][1] = -q * ny + r * nx;
        double hR = h_riemann_int[edge_idx][1];
        double cR = std::sqrt(gravity * hR);
        lambda_plus_int[edge_idx][1] = lambdaPlus(cR);
        lambda_minus_int[edge_idx][1] = lambdaMinus(cR);
		// === Right-hand side
		lambda_max_int[edge_idx] = std::max(lambdaPlus(cL) , lambdaPlus(cR));
        
	}
    //==== internal edges
    return EXIT_SUCCESS;
}
