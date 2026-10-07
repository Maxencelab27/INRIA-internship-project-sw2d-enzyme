
#include <cmath>

#include "swLiswes.h"
#include "../swParent.h"
#include "../common/swHyp.h"
#include "../../../common/param.h"
#include "../../../common/swLog.h"

void swLiswes::computeMusclSlopes(const std::vector<double>& cell_x, const std::vector<double>& cell_y, const std::vector<std::vector<size_t>>& cell_neigh_cell_idx,
                               const std::vector<double>& edge_int_midpoint_x, const std::vector<double>& edge_int_midpoint_y,
                               const std::vector<size_t>& edge_int_celL_idx, const std::vector<size_t>& edge_int_celR_idx,
                               const std::vector<double>& cell_averaged_field, std::vector<double>& slope_x, std::vector<double>& slope_y)
/**
 * @brief Computes the MUSCL coefficients in the x- and y-directions

 * @param[in] cell_x[i]: x-coordinate of cell with index i
 * @param[in] cell_y[i]: y-coordinate of cell with index i
 * @param[in] cell_neigh_cell_idx[i][j] : idex of the j-th neighbour cell of the i-th cell
 * @param[in] edge_int_midpoint_x[i]: x-coordinate of the i-th internal edge
 * @param[in] edge_int_midpoint_y[i]: y-coordinate of the i-th internal edge
 * @param[in] edge_int_celL_idx[i]: index othe LHS cell for the i-th internal edge
 * @param[in] edge_int_celR_idx[i]: index othe RHS cell for the i-th internal edge
 * @param[in] cell_averaged_field[i] : the average value in cell idx i of the variable to be reconstructed over the computational cells
 * @param[in] cell_neigh_cell_idx[i][j] : idex of the j-th neighbour cell of the i-th cell
 * @param[out] slope_x[i]: slope in the x-direction for cell index i
 * @param[out] slope_y[i]: slope in the y-direction for thell with index i
 */
{
    //-- Step#1: compute 2nd-order-like gradients
    size_t cell_nb = cell_averaged_field.size();
    for (size_t idx_cell = 0;idx_cell < cell_nb; idx_cell++)
    {
        slope_x[idx_cell] = 0.0;
        slope_y[idx_cell] = 0.0;
        double weight_sum = 0.0;
        if (cell_neigh_cell_idx[idx_cell].size() > 1)
        {
            //---- Compute individual contributions of neighbour cells - only if more than 1 neighbour
            for (size_t cell_iterator = 0; cell_iterator < cell_neigh_cell_idx[idx_cell].size(); cell_iterator++)
            {
                size_t j = cell_neigh_cell_idx[idx_cell][cell_iterator];
                double dx = cell_x[j] - cell_x[idx_cell];
                double dy = cell_y[j] - cell_y[idx_cell];
                double dist = std::sqrt(dx * dx + dy * dy);
                double weight = 1.0 / dist;
                slope_x[idx_cell] += (cell_averaged_field[j] - cell_averaged_field[idx_cell]) * dx / dist / dist * weight;
                slope_y[idx_cell] += (cell_averaged_field[j] - cell_averaged_field[idx_cell]) * dy / dist / dist * weight;
                weight_sum += weight;
            }
            //==== Compute individual contributions of neighbour cells - only if more than 1 neighbour
            //---- Normalize
            slope_x[idx_cell] /= weight_sum;
            slope_y[idx_cell] /= weight_sum;
            //==== Normalize
        }
    }
    //== Step#1: compute 2nd-order-like gradients

    //-- Step#2: apply limiting
    for (size_t idx_edge = 0; idx_edge < edge_int_midpoint_x.size(); idx_edge++)
    {
        //---- Compute the min & max of field on each side
        size_t idx_celL = edge_int_celL_idx[idx_edge];
        size_t idx_celR = edge_int_celR_idx[idx_edge];
        double field_max = std::max(cell_averaged_field[idx_celL], cell_averaged_field[idx_celR]);
        double field_min = std::min(cell_averaged_field[idx_celL], cell_averaged_field[idx_celR]);
        //==== Compute the min & max of field on each side

        //---- LHS: compute reconstructed value and apply limiting if necessary
        double field_reconst = cell_averaged_field[idx_celL] + (edge_int_midpoint_x[idx_edge] - cell_x[idx_celL]) * slope_x[idx_celL]
                                                             + (edge_int_midpoint_y[idx_edge] - cell_y[idx_celL]) * slope_y[idx_celL];
        if (field_reconst < field_min) // Undershooting
        {
            double delta_f_reconst = field_reconst - cell_averaged_field[idx_celL];
            double delta_f_target  = field_min - cell_averaged_field[idx_celL];
            double slope_factor = delta_f_target / delta_f_reconst;
            slope_x[idx_celL] *= slope_factor;
            slope_y[idx_celL] *= slope_factor;
        }
        if (field_reconst > field_max) // Overshooting
        {
            double delta_f_reconst = field_reconst - cell_averaged_field[idx_celL];
            double delta_f_target = field_max - cell_averaged_field[idx_celL];
            double slope_factor = delta_f_target / delta_f_reconst;
            slope_x[idx_celL] *= slope_factor;
            slope_y[idx_celL] *= slope_factor;
        }
        //==== LHS: compute reconstructed value and apply limiting if necessary

        //---- RHS: compute reconstructed value and apply limiting if necessary
        field_reconst = cell_averaged_field[idx_celR] + (edge_int_midpoint_x[idx_edge] - cell_x[idx_celR]) * slope_x[idx_celR]
                                                      + (edge_int_midpoint_y[idx_edge] - cell_y[idx_celR]) * slope_y[idx_celR];
        if (field_reconst < field_min) // Undershooting
        {
            double delta_f_reconst = field_reconst - cell_averaged_field[idx_celR];
            double delta_f_target = field_min - cell_averaged_field[idx_celR];
            double slope_factor = delta_f_target / delta_f_reconst;
            slope_x[idx_celR] *= slope_factor;
            slope_y[idx_celR] *= slope_factor;
        }
        if (field_reconst > field_max) // Overshooting
        {
            double delta_f_reconst = field_reconst - cell_averaged_field[idx_celR];
            double delta_f_target = field_max - cell_averaged_field[idx_celR];
            double slope_factor = delta_f_target / delta_f_reconst;
            slope_x[idx_celR] *= slope_factor;
            slope_y[idx_celR] *= slope_factor;
        }
        //==== RHS: compute reconstructed value and apply limiting if necessary


    }
}


bool swLiswes::hypMusclEvrComputeLRStates(void)
/**
 * @brief Computes the left and right states ( z, q, r, lambda+/- ) on every edge.
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 *   - all interface-based Riemann states
 */
{
    //-- Compute free surface slopes
    std::vector<double> z_cell = {};
    for (size_t cell_idx = 0; cell_idx < cell_nb; cell_idx++)
    {
        z_cell.push_back(water_depth[cell_idx] + cell_z_b[cell_idx]);
    }
    std::vector<double> dz_dx = {};
    dz_dx.resize(cell_nb);
    std::vector<double> dz_dy = {};
    dz_dy.resize(cell_nb);
    computeMusclSlopes(cell_centre_x, cell_centre_y, cell_neigh_cell_idx,
                              edge_int_midpoint_x, edge_int_midpoint_y, 
                              edge_int_celL_idx, edge_int_celR_idx,
                              z_cell, dz_dx, dz_dy);
    //==Compute free surface slopes
    
    //---- Boundary edges
    for (size_t edge_idx = 0 ; edge_idx < edge_bnd_nb ; edge_idx++)
    {
        size_t idx_cell = edge_bnd_cell_idx[edge_idx];
        double nx = edge_bnd_norm_x[edge_idx];
        double ny = edge_bnd_norm_y[edge_idx];
        double q = unit_discharge_x[idx_cell];
        double r = unit_discharge_y[idx_cell];
        h_riemann_bnd[edge_idx] = water_depth[idx_cell];
        z_riemann_bnd[edge_idx] = water_depth[idx_cell] + cell_z_b[idx_cell];
        q_riemann_bnd[edge_idx] = q * nx + r * ny;
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
        //-- Prepare variables needed on both sides of interface
        double nx = edge_int_norm_x[edge_idx];
        double ny = edge_int_norm_y[edge_idx];
        size_t idx_celL = edge_int_celL_idx[edge_idx];
        size_t idx_celR = edge_int_celR_idx[edge_idx];
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
        double u_local = q_local / std::max(small_real, h_local);
        double v_local = r_local / std::max(small_real, h_local);
        //------ Lambda 1 (u-c)
        double lambda_1 = std::max(0.0, lambdaMinus(c_local));
        double delta_x_1_L =  edge_int_midpoint_x[edge_idx] - cell_centre_x[idx_celL]
                            - 0.5 * lambda_1 * dt * nx; // delta_x is the x-distance of the middle of the depdomain to the cell centre
        double delta_y_1_L =  edge_int_midpoint_y[edge_idx] - cell_centre_y[idx_celL]
                            - 0.5 * lambda_1 * dt * ny; // delta_y is the y-distance of the middle of the depdomain to the cell centre
        double z_evr_1_L = z_cell[idx_celL] + delta_x_1_L * dz_dx[idx_celL] + delta_y_1_L * dz_dy[idx_celL];
        double h_evr_1_L = std::max(0.0, z_evr_1_L - z_b_L);
        //------ Lambda 1 (u-c)
        //------ Lambda 2 (u+c)
        double lambda_2 = std::max(0.0, lambdaPlus(c_local));
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
        double cL = std::sqrt(gravity * hL);
        lambda_plus_int[edge_idx][0] = lambdaPlus(cL);
        lambda_minus_int[edge_idx][0] = lambdaMinus(cL);
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
        u_local = q_local / std::max(small_real, h_local);
        v_local = r_local / std::max(small_real, h_local);
        //------ Lambda 1 (u-c)
        lambda_1 = std::min(0.0, lambdaMinus(c_local));
        double delta_x_1_R =  edge_int_midpoint_x[edge_idx] - cell_centre_x[idx_celR]
                             - 0.5 * lambda_1 * dt * nx; // delta_x is the x-distance of the middle of the depdomain to the cell centre
        double delta_y_1_R =  edge_int_midpoint_y[edge_idx] - cell_centre_y[idx_celR]
                             - 0.5 * lambda_1 * dt * ny; // delta_y is the y-distance of the middle of the depdomain to the cell centre
        double z_evr_1_R = z_cell[idx_celR] + delta_x_1_R * dz_dx[idx_celR] + delta_y_1_R * dz_dy[idx_celR];
        double h_evr_1_R = std::max(0.0, z_evr_1_R - z_b_R);
        //------ Lambda 1 (u-c)
        //------ Lambda 2 (u+c)
        lambda_2 = std::min(0.0, lambdaPlus(c_local));
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
        double cR = std::sqrt(gravity * hR);
        lambda_plus_int[edge_idx][1] = lambdaPlus(cR);
        lambda_minus_int[edge_idx][1] = lambdaMinus(cR);
		//== Right-hand side
		lambda_max_int[edge_idx] = std::max(lambdaPlus(cL) , lambdaPlus(cR));
}
    //==== internal edges
    return EXIT_SUCCESS;
}
