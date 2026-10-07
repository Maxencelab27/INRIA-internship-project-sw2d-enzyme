/**
 * @file swSpHypMusclEvr_E.cpp
 * @brief Differentiated (Enzyme) counterpart of swSpHypMusclEvr.cpp.
 *
 * Added for automatic differentiation; not part of native SW2D.
 * Declarations live in swSp.h (AD / ENZYME block).
 */

#include <cmath>
#include <algorithm>

#include "swSpE.h"
#include "../swParent.h"
#include "../common/swHyp.h"
#include "../../../common/param.h"
#include "../../../common/swLog.h"
#include "swSpDesingVelocity.h"

using swSpDesing::desingVelocity;
using swSpDesing::kHDesing;

bool swSpE::hypMusclEvrComputeLRStates_E(void)
/**
 * @brief Computes the left and right states ( z, q, r, lambda+/- ) on every edge.
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 *   - all interface-based Riemann states
 */
{
    //-- Compute free surface slopes

    

   
    for (size_t cell_idx = 0; cell_idx < cell_nb; cell_idx++)
    {
        this->z_cell[cell_idx] = water_depth[cell_idx] + cell_z_b[cell_idx];
    }


    computeMusclSlopes(cell_centre_x, cell_centre_y, cell_neigh_cell_idx,
                       edge_int_midpoint_x, edge_int_midpoint_y, 
                       edge_int_celL_idx, edge_int_celR_idx,
                       this->z_cell, this->dz_dx, this->dz_dy);

   
    
    //==Compute free surface slopes
    
    //---- Boundary edges
    for (size_t edge_idx = 0 ; edge_idx < edge_bnd_nb ; edge_idx++)
    {
        double nx = edge_bnd_norm_x[edge_idx];
        double ny = edge_bnd_norm_y[edge_idx];
        size_t idx_cell = edge_bnd_cell_idx[edge_idx];
        double bR = boussinesq[idx_cell];
        double q = unit_discharge_x[idx_cell];
        double r = unit_discharge_y[idx_cell];
        this->h_riemann_bnd[edge_idx] = water_depth[idx_cell];
        this->z_riemann_bnd[edge_idx] = water_depth[idx_cell] + cell_z_b[idx_cell];
        this->q_riemann_bnd[edge_idx] = q * nx + r * ny;
        this->r_riemann_bnd[edge_idx] = -q * ny + r * nx;
        double hR = h_riemann_bnd[edge_idx];
        double uR = desingVelocity(q_riemann_bnd[edge_idx], hR, kHDesing);
        double cR = std::sqrt(gravity * hR);
        lambda_minus_bnd[edge_idx] = lambdaMinus(bR, cR, uR);
        lambda_plus_bnd[edge_idx] = lambdaPlus(bR, cR, uR);
		lambda_max_bnd[edge_idx] = std::max(std::abs(lambda_minus_bnd[edge_idx]), std::abs(lambda_plus_bnd[edge_idx]));
	}
    //==== Boundary edges
	
    //---- Internal edges
    for (size_t edge_idx = 0 ; edge_idx < edge_int_nb ; edge_idx++)
    {
        //-- Prepare variables needed on both sides of interface
        double nx = edge_int_norm_x[edge_idx];
        double ny = edge_int_norm_y[edge_idx];
        size_t idx_celL = edge_int_celL_idx[edge_idx];
        size_t idx_celR = edge_int_celR_idx[edge_idx];
        double bL = boussinesq[idx_celL];
        double bR = boussinesq[idx_celR];
        double z_b_L = cell_z_b[idx_celL];
        double z_b_R = cell_z_b[idx_celR];
        //== Prepare variables needed on both sides of interface
        
		//-- Left-hand side
        //---- Convert (q,r) & surface slopes to local coordinate system
		double q = unit_discharge_x[idx_celL];
		double r = unit_discharge_y[idx_celL];
        double h_local = water_depth[idx_celL];
        double q_local = q * nx + r * ny;
        double r_local = r * nx - q * ny;
        //====  Convert (q,r) & surface slopes to local coordinate system

        //---- Compute lambda(1/2)*dt and infer average values
        double c_local = std::sqrt(gravity * std::max(small_real, h_local));
        double u_local = desingVelocity(q_local, h_local, kHDesing);
        double v_local = desingVelocity(r_local, h_local, kHDesing);
        //------ Lambda 1 (u-c)
        double lambda_1 = std::max(0.0, lambdaMinus(bL, c_local, u_local));
        double delta_x_1_L =  edge_int_midpoint_x[edge_idx] - cell_centre_x[idx_celL]
                            - 0.5 * lambda_1 * dt * nx; // delta_x is the x-distance of the middle of the depdomain to the cell centre
        double delta_y_1_L =  edge_int_midpoint_y[edge_idx] - cell_centre_y[idx_celL]
                            - 0.5 * lambda_1 * dt * ny; // delta_y is the y-distance of the middle of the depdomain to the cell centre
        double z_evr_1_L = z_cell[idx_celL] + delta_x_1_L * dz_dx[idx_celL] + delta_y_1_L * dz_dy[idx_celL];
        double h_evr_1_L = std::max(0.0, z_evr_1_L - z_b_L);
        //------ Lambda 1 (u-c)
        //------ Lambda 2 (u+c)
        double lambda_2 = std::max(0.0, lambdaPlus(bL, c_local, u_local));
        double delta_x_2_L =  edge_int_midpoint_x[edge_idx] - cell_centre_x[idx_celL]
                             - 0.5 * lambda_2 * dt * nx; // delta_x is the x-distance of the middle of the depdomain to the cell centre
        double delta_y_2_L =  edge_int_midpoint_y[edge_idx] - cell_centre_y[idx_celL]
                            - 0.5 * lambda_2 * dt * ny; // delta_y is the y-distance of the middle of the depdomain to the cell centre
        double z_evr_2_L = z_cell[idx_celL] + delta_x_2_L * dz_dx[idx_celL] + delta_y_2_L * dz_dy[idx_celL];
        double h_evr_2_L = std::max(0.0, z_evr_2_L - z_b_L);
        //------ Lambda 1 (u+c)
        //== Compute lambda(1/2)*dt and infer average values

        h_riemann_int[edge_idx][0] = 0.5 * (h_evr_1_L + h_evr_2_L);
        z_riemann_int[edge_idx][0] = 0.5 * (h_evr_1_L + h_evr_2_L) + cell_z_b[idx_celL];
        q_riemann_int[edge_idx][0] = 0.5 * (h_evr_1_L + h_evr_2_L) * u_local;
        r_riemann_int[edge_idx][0] = 0.5 * (h_evr_1_L + h_evr_2_L) * v_local;
        double hL = h_riemann_int[edge_idx][0];
        double uL = desingVelocity(q_riemann_int[edge_idx][0], hL, kHDesing);
        double cL = std::sqrt(gravity * hL);
        lambda_minus_int[edge_idx][0] = lambdaMinus(bL, cL, uL);
        lambda_plus_int[edge_idx][0] = lambdaPlus(bL, cL, uL);
		//== Left-hand side
		



		//-- Right-hand side
        //---- Convert (q,r) & surface slopes to local coordinate system
        q = unit_discharge_x[idx_celR];
        r = unit_discharge_y[idx_celR];
        h_local = water_depth[idx_celR];
        q_local = unit_discharge_x[idx_celR] * nx + unit_discharge_y[idx_celR] * ny;
        r_local = unit_discharge_y[idx_celR] * nx - unit_discharge_x[idx_celR] * ny;
        //====  Convert (q,r) & surface slopes to local coordinate system

        //---- Compute lambda(1/2)*dt and infer average values
        c_local = std::sqrt(gravity * std::max(small_real, h_local));
        u_local = desingVelocity(q_local, h_local, kHDesing);
        v_local = desingVelocity(r_local, h_local, kHDesing);
        //------ Lambda 1 (u-c)
        lambda_1 = std::min(0.0, lambdaMinus(bR, c_local, u_local));
        double delta_x_1_R =  edge_int_midpoint_x[edge_idx] - cell_centre_x[idx_celR]
                             - 0.5 * lambda_1 * dt * nx; // delta_x is the x-distance of the middle of the depdomain to the cell centre
        double delta_y_1_R =  edge_int_midpoint_y[edge_idx] - cell_centre_y[idx_celR]
                             - 0.5 * lambda_1 * dt * ny; // delta_y is the y-distance of the middle of the depdomain to the cell centre
        double z_evr_1_R = z_cell[idx_celR] + delta_x_1_R * dz_dx[idx_celR] + delta_y_1_R * dz_dy[idx_celR];
        double h_evr_1_R = std::max(0.0, z_evr_1_R - z_b_R);
        //------ Lambda 1 (u-c)
        //------ Lambda 2 (u+c)
        lambda_2 = std::min(0.0, lambdaPlus(bR, c_local, u_local));
        double delta_x_2_R =  edge_int_midpoint_x[edge_idx] - cell_centre_x[idx_celR]
                             - 0.5 * lambda_2 * dt * nx; // delta_x is the x-distance of the middle of the depdomain to the cell centre
        double delta_y_2_R =  edge_int_midpoint_y[edge_idx] - cell_centre_y[idx_celR]
                            - 0.5 * lambda_2 * dt * ny; // delta_y is the y-distance of the middle of the depdomain to the cell centre
        double z_evr_2_R = z_cell[idx_celR] + delta_x_2_R * dz_dx[idx_celR] + delta_y_2_R * dz_dy[idx_celR];
        double h_evr_2_R = std::max(0.0, z_evr_2_R - z_b_R);
        //------ Lambda 1 (u+c)
        //== Compute lambda(1/2)*dt and infer average values

        h_riemann_int[edge_idx][1] = 0.5 * (h_evr_1_R + h_evr_2_R);
        z_riemann_int[edge_idx][1] = 0.5 * (h_evr_1_R + h_evr_2_R) + cell_z_b[idx_celR];
        q_riemann_int[edge_idx][1] = 0.5 * (h_evr_1_R + h_evr_2_R) * u_local;
        r_riemann_int[edge_idx][1] = 0.5 * (h_evr_1_R + h_evr_2_R) * v_local;
        double hR = h_riemann_int[edge_idx][1];
        double uR = desingVelocity(q_riemann_int[edge_idx][1], hR, kHDesing);
        double cR = std::sqrt(gravity * hR);
        lambda_minus_int[edge_idx][1] = lambdaMinus(bR, cR, uR);
        lambda_plus_int[edge_idx][1] = lambdaPlus(bR, cR, uR);
		//== Right-hand side
        // Use { } in max hereafter to use more thn 2 variables
        double m1 = std::max(std::abs(lambda_minus_int[edge_idx][0]), std::abs(lambda_plus_int[edge_idx][0]));
        double m2 = std::max(std::abs(lambda_minus_int[edge_idx][1]), std::abs(lambda_plus_int[edge_idx][1]));
        lambda_max_int[edge_idx] = std::max(m1, m2);
		/*lambda_max_int[edge_idx] = std::max({std::abs(lambda_minus_int[edge_idx][0]),
                                             std::abs(lambda_plus_int[edge_idx][0]), 
                                             std::abs(lambda_minus_int[edge_idx][1]), 
                                             std::abs(lambda_plus_int[edge_idx][1])});*/
}
    //==== internal edges
    return EXIT_SUCCESS;
}
