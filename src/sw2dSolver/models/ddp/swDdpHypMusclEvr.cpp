
#include <cmath>

#include "swDdp.h"
#include "../swParent.h"
#include "../common/swHyp.h"
#include "../../../common/param.h"
#include "../../../common/swLog.h"


/** @todo: adapt code */
bool swDdp::computeMusclSlopes(const std::vector<double>& cell_x, const std::vector<double>& cell_y, const std::vector<std::vector<size_t>>& cell_neigh_cell_idx,
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
    for (size_t idx_edge_int= 0; idx_edge_int< edge_int_midpoint_x.size(); idx_edge_int++)
    {
        //---- Compute the min & max of field on each side
        size_t idx_celL = edge_int_celL_idx[idx_edge_int];
        size_t idx_celR = edge_int_celR_idx[idx_edge_int];
        double field_max = std::max(cell_averaged_field[idx_celL], cell_averaged_field[idx_celR]);
        double field_min = std::min(cell_averaged_field[idx_celL], cell_averaged_field[idx_celR]);
        //==== Compute the min & max of field on each side

        //---- LHS: compute reconstructed value and apply limiting if necessary
        double field_reconst = cell_averaged_field[idx_celL] + (edge_int_midpoint_x[idx_edge_int] - cell_x[idx_celL]) * slope_x[idx_celL]
                                                             + (edge_int_midpoint_y[idx_edge_int] - cell_y[idx_celL]) * slope_y[idx_celL];
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
        field_reconst = cell_averaged_field[idx_celR] + (edge_int_midpoint_x[idx_edge_int] - cell_x[idx_celR]) * slope_x[idx_celR]
                                                      + (edge_int_midpoint_y[idx_edge_int] - cell_y[idx_celR]) * slope_y[idx_celR];
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
    return EXIT_SUCCESS;
}


bool swDdp::hypMusclEvrComputeLRStates(void)
/**
 * @brief Computes the left and right states ( z, q, r, lambda+/- ) on every edge.
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 *   - all interface-based Riemann states
 */
{
    return EXIT_SUCCESS;
}
