#include <cmath>

#include "swDdp.h"
#include "../../../common/param.h"
#include "../../../common/swLog.h"


bool swDdp::hypGodunovComputeLRStates(void)
/**
 * @brief Computes the left and right states ( z, q, r, lambda+/- ) on every edge
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 *   - all interface-based Riemann states
 */
{
    ///-- Boundary edges
	for (size_t idx_edge_bnd = 0; idx_edge_bnd < edge_bnd_nb; idx_edge_bnd++)
    {
        double normX_edge = edge_bnd_norm_x[idx_edge_bnd];
        double normY_edge = edge_bnd_norm_y[idx_edge_bnd];

        //-- get parameters in the computational cell
        size_t idx_cell = edge_bnd_cell_idx[idx_edge_bnd];
        double theta_in = thetaW[idx_cell];
        double z_in = zOfThetaW(idx_cell, theta_in);
        double h_in = std::max(z_in - cell_z_tab[idx_cell][0], 0.);
        theta_riemann_bnd[idx_edge_bnd] = h_in;
        //== get parameters in the computational cell     

        //-- recompute hydrodynamic values in the computational cell
        double theta_u = thetaW_velocity_x[idx_cell];
        double theta_v = thetaW_velocity_y[idx_cell];
        theta_u_riemann_bnd[idx_edge_bnd] = theta_u * normX_edge + theta_v * normY_edge;
        theta_v_riemann_bnd[idx_edge_bnd] = -theta_u * normY_edge + theta_v * normX_edge;
        if (h_in < h_min)
        {
            theta_u_riemann_bnd[idx_edge_bnd] = 0.;
            theta_v_riemann_bnd[idx_edge_bnd] = 0.;
        }
        z_riemann_bnd[idx_edge_bnd] = z_in;
        //== recompute hydrodynamic values in the computational cell

        //-- compute characteristic variables
        double phi_omega = phiWOfZ(idx_cell, z_in);
        /// On boundaries, phi_gamma is forced to be == phi_omega (and so are thetas)
        double x_phi = 1.0; // phi_omega / std::max(phi_gamma, small_real);
        double x_theta = 1.0; // theta_omega / std::max(theta_gamma, small_real);

        double c_in = pressureWaveCelerity(theta_in, phi_omega);
        double b_in = boussinesq[idx_cell];
        double u_in = theta_u_riemann_bnd[idx_edge_bnd] / std::max(theta_in, small_real);
        lambda_plus_bnd[idx_edge_bnd] = lambdaPlus(u_in, c_in, x_phi, x_theta, b_in);
        lambda_minus_bnd[idx_edge_bnd] = lambdaMinus(u_in, c_in, x_phi, x_theta, b_in);
        lambda_max_bnd[idx_edge_bnd] = std::max(std::abs(lambda_plus_bnd[idx_edge_bnd]),std::abs(lambda_minus_bnd[idx_edge_bnd]));
    }
    ///== Boundary edges

    ///-- Internal edges
    for (size_t idx_edge_int = 0; idx_edge_int < edge_int_nb; idx_edge_int++)
    {
        double normX_edge = edge_int_norm_x[idx_edge_int];
        double normY_edge = edge_int_norm_y[idx_edge_int];

        //-- Left-hand side
        //---- get parameters in the computational cell
        size_t idx_celL = edge_int_celL_idx[idx_edge_int];
        double theta_L = thetaW[idx_celL];
        double z_L = zOfThetaW(idx_celL, theta_L);
        double h_L = std::max(z_L - cell_z_tab[idx_celL][0], 0.);
        theta_riemann_int[idx_edge_int][0] = theta_L;
        //== get parameters in the computational cell     

        //-- recompute hydrodynamic values in the computational cell
        double theta_u_L = thetaW_velocity_x[idx_celL];
        double theta_v_L = thetaW_velocity_y[idx_celL];
        theta_u_riemann_int[idx_edge_int][0] = theta_u_L * normX_edge + theta_v_L * normY_edge;
        theta_v_riemann_int[idx_edge_int][0] = -theta_u_L * normY_edge + theta_v_L * normX_edge;
        z_riemann_int[idx_edge_int][0] = z_L;
        //== recompute hydrodynamic values in the computational cell

        //-- compute characteristic variables
        double phi_omega = phiWOfZ(idx_celL, z_L);
        double phi_gamma_L = phiGIntOfZ(idx_edge_int, z_L);
        double theta_omega = thetaWOfZ(idx_celL, z_L);
        double theta_gamma_L = thetaGIntOfZ(idx_edge_int, z_L);
        double x_phi = 1.0;
        double x_theta = 1.0;
        if (theta_gamma_L * theta_omega > small_real) // avoid empty interface
        {
            if (h_L < h_min)
            {
                theta_u_riemann_int[idx_edge_int][0] = 0.;
                theta_v_riemann_int[idx_edge_int][0] = 0.;
                double z_min = cell_z_tab[idx_celL][0] + h_min;
                phi_omega = phiWOfZ(idx_celL, z_min);
                phi_gamma_L = phiGIntOfZ(idx_edge_int, z_min);
                theta_omega = thetaWOfZ(idx_celL, z_min);
                theta_gamma_L = thetaGIntOfZ(idx_edge_int, z_min);
            }
            x_phi = phi_omega / std::max(phi_gamma_L, small_real);
            x_theta = theta_omega / std::max(theta_gamma_L, small_real); 
        }
        double u_L = theta_u_riemann_int[idx_edge_int][0] / std::max(theta_L, small_real);
        double c_L = pressureWaveCelerity(theta_L, phi_omega);
        double b_L = boussinesq[idx_celL];
        lambda_plus_int[idx_edge_int][0] = lambdaPlus(u_L, c_L, x_phi, x_theta, b_L);
        lambda_minus_int[idx_edge_int][0] = lambdaMinus(u_L, c_L, x_phi, x_theta, b_L);
        //== Left-hand side

        //-- Right-hand side
        //---- get parameters in the computational cell
        size_t idx_celR = edge_int_celR_idx[idx_edge_int];
        double theta_R = thetaW[idx_celR];
        double z_R = zOfThetaW(idx_celR, theta_R);
        double h_R = std::max(z_R - cell_z_tab[idx_celR][0], 0.);
        theta_riemann_int[idx_edge_int][1] = theta_R;
        //== get parameters in the computational cell     

        //-- recompute hydrodynamic values in the computational cell
        double theta_u_R = thetaW_velocity_x[idx_celR];
        double theta_v_R = thetaW_velocity_y[idx_celR];
        theta_u_riemann_int[idx_edge_int][1] = theta_u_R * normX_edge + theta_v_R * normY_edge;
        theta_v_riemann_int[idx_edge_int][1] = -theta_u_R * normY_edge + theta_v_R * normX_edge;
        z_riemann_int[idx_edge_int][1] = z_R;
        //== recompute hydrodynamic values in the computational cell
        //-- compute characteristic variables
        phi_omega = phiWOfZ(idx_celR, z_R);
        double phi_gamma_R = phiGIntOfZ(idx_edge_int, z_R);
        theta_omega = thetaWOfZ(idx_celR, z_R);
        double theta_gamma_R = thetaGIntOfZ(idx_edge_int, z_R);
        x_phi = 1.0;
        x_theta = 1.0;
        if (theta_gamma_R * theta_omega > small_real) // avoid empty interface
        {
            if (h_R < h_min)
            {
                theta_u_riemann_int[idx_edge_int][1] = 0.;
                theta_v_riemann_int[idx_edge_int][1] = 0.;
                double z_min = cell_z_tab[idx_celR][1] + h_min;
                phi_omega = phiWOfZ(idx_celR, z_min);
                phi_gamma_R = phiGIntOfZ(idx_edge_int, z_min);
                theta_omega = thetaWOfZ(idx_celR, z_min);
                theta_gamma_R = thetaGIntOfZ(idx_edge_int, z_min);
            }
            x_phi = phi_omega / std::max(phi_gamma_R, small_real);
            x_theta = theta_omega / std::max(theta_gamma_R, small_real); 
        }
        double u_R = theta_u_riemann_int[idx_edge_int][1] / std::max(theta_R, small_real);
        double c_R = pressureWaveCelerity(theta_R, phi_omega);
        double b_R = boussinesq[idx_celR];
        lambda_plus_int[idx_edge_int][1] = lambdaPlus(u_R, c_R, x_phi, x_theta, b_R);
        lambda_minus_int[idx_edge_int][1] = lambdaMinus(u_R, c_R, x_phi, x_theta, b_R);
        //== Right-hand side
        
        // Use { } in max hereafter to use more than 2 variables
		lambda_max_int[idx_edge_int] = std::max({std::abs(lambda_minus_int[idx_edge_int][0]),
                                             std::abs(lambda_plus_int[idx_edge_int][0]), 
                                             std::abs(lambda_minus_int[idx_edge_int][1]), 
                                             std::abs(lambda_plus_int[idx_edge_int][1])});
    }
    ///== Internal edges
    return EXIT_SUCCESS;
}