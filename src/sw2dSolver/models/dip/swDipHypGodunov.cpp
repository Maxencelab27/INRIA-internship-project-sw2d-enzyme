#include <cmath>

#include "swDip.h"
#include "../../../common/param.h"
#include "../../../common/swLog.h"


bool swDip::hypGodunovComputeLRStates(void)
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
        double bR = boussinesq[id_cell];
        double phiG = phiG_bnd[edge_idx];
        double phiR = phiW[id_cell];
		double q = unit_discharge_x[id_cell];
		double r = unit_discharge_y[id_cell];
        h_riemann_bnd[edge_idx] = water_depth[id_cell];
        z_riemann_bnd[edge_idx] = water_depth[id_cell] + cell_z_b[id_cell];
        q_riemann_bnd[edge_idx] = q * nx  + r * ny;
        r_riemann_bnd[edge_idx] = -q * ny + r * nx;
        double hR = h_riemann_bnd[edge_idx];
        double uR = q_riemann_bnd[edge_idx] / std::max(hR, small_real);
        double cR = std::sqrt(gravity * hR);
        lambda_minus_bnd[edge_idx] = lambdaMinus(bR, phiG, phiR, cR, uR);
        lambda_plus_bnd[edge_idx] = lambdaPlus(bR, phiG, phiR, cR, uR);
		lambda_max_bnd[edge_idx] = std::max(std::abs(lambda_plus_bnd[edge_idx]),std::abs(lambda_minus_bnd[edge_idx]));
	}
    //==== Boundary edges
	
    //---- Internal edges
    for (size_t edge_idx = 0 ; edge_idx < edge_int_nb ; edge_idx++)
    {
        double nx = edge_int_norm_x[edge_idx];
        double ny = edge_int_norm_y[edge_idx];
        double phiG = phiG_int[edge_idx];
        
		// --- Left-hand side
        size_t id_celL = edge_int_celL_idx[edge_idx];
        double bL = boussinesq[id_celL];
        double phiL = phiW[id_celL];
		double q = unit_discharge_x[id_celL];
		double r = unit_discharge_y[id_celL];
        h_riemann_int[edge_idx][0] = water_depth[id_celL];
        z_riemann_int[edge_idx][0] = water_depth[id_celL] + cell_z_b[id_celL];
        q_riemann_int[edge_idx][0] = q * nx  + r * ny;
        r_riemann_int[edge_idx][0] = -q * ny + r * nx;
        double hL = h_riemann_int[edge_idx][0];
        double uL = q_riemann_int[edge_idx][0] / std::max(hL, small_real);
        double cL = std::sqrt(gravity * hL);
        lambda_minus_int[edge_idx][0] = lambdaMinus(bL, phiG, phiL, cL, uL);
        lambda_plus_int[edge_idx][0] = lambdaPlus(bL, phiG, phiL, cL, uL);
		// === Left-hand side
		
		// --- Right-hand side
        size_t id_celR = edge_int_celR_idx[edge_idx];
        double bR = boussinesq[id_celR];
        double phiR = phiW[id_celR];
		q = unit_discharge_x[id_celR];
		r = unit_discharge_y[id_celR];
        h_riemann_int[edge_idx][1] = water_depth[id_celR];
        z_riemann_int[edge_idx][1] = water_depth[id_celR] + cell_z_b[id_celR];
        q_riemann_int[edge_idx][1] = q * nx  + r * ny;
        r_riemann_int[edge_idx][1] = -q * ny + r * nx;
        double hR = h_riemann_int[edge_idx][1];
        double uR = q_riemann_int[edge_idx][1] / std::max(hR, small_real);
        double cR = std::sqrt(gravity * hR);
        lambda_minus_int[edge_idx][1] = lambdaMinus(bR, phiG, phiR, cR, uR);
        lambda_plus_int[edge_idx][1] = lambdaPlus(bR, phiG, phiR, cR, uR);
		// === Right-hand side
        // Use { } in max hereafter to use more than 2 variables
		lambda_max_int[edge_idx] = std::max({std::abs(lambda_minus_int[edge_idx][0]),
                                             std::abs(lambda_plus_int[edge_idx][0]), 
                                             std::abs(lambda_minus_int[edge_idx][1]), 
                                             std::abs(lambda_plus_int[edge_idx][1])});
    }
    //==== internal edges
    return EXIT_SUCCESS;
}
