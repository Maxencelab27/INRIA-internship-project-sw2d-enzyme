#include <algorithm>
#include <cmath>
#include <vector>

#include "./swDdp.h"
#include "../common/swSource.h"
#include "../common/swBoundaryConditionTimeSeries.h"
#include "../../../common/swBoundaryCondition.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"
#include "../../../common/param.h"


bool swDdp::hypBoundaryFlux(const size_t& idx_edge_bnd, const double& bR, const double& zbR, const double& zR, 
                           const double& uR, const double& vR, const double& lambda_plus_R, 
                           std::string const& bc_type, const double& bc_val,
                           double& flux1, double& flux2, double& flux3, double& source_right)
/**
 * @brief Computes fluxes at boundary edges for the DDP models
 * @param[in] idx_edge_bnd  : boundary edge index
 * @param[in] bR  : Boussinesq coefficient on the RHS of the interface
 * @param[in] zbR : bottom level in the cell on the RHS of the interface
 * @param[in] zR  : free surface elevation in the cell on the RHS of the interface
 * @param[in] uR  : normal velocity in the cell on the RHS of the interface
 * @param[in] vR  : tangent velocity in the cell on the RHS of the interface
 * @param[in] lambda_plus_R : lambda+ value in the cell on the RHS of the interface
 * @param[in] bc_type : type of boundary condition on the boundary interface ('h', 'q', 'fr', etc.)
 * @param[in] bc_val : value for the boundary condition of type bc_type
 * @param[out] flux1 : mass flux in the direction normal to the interface
 * @param[out] flux2 : momentum flux in the direction normal to the interface (specific force per unit width)
 * @param[out] flux3 : momentum flux in the direction tangent to the interface
 * @param[out] source_right: geometric source term on the RHS of the interface
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 */
    /* Formulae are designed for a domain on the right-hand side of the boundary */
{
    //-- compute characteristic variables
    size_t idx_cell = edge_bnd_cell_idx[idx_edge_bnd];
    double pressure_force = pressureForceOmegaOfZ(idx_cell, zR);    
    double phi_omega = phiWOfZ(idx_cell, zR);
    double phi_gamma = phi_omega; // On boundaries, phi_gamma is forced to be == phi_omega
    double x_phi = phi_omega / std::max(phi_gamma, small_real);
 

    double theta_omega = thetaWOfZ(idx_cell, zR);
    double theta_gamma = theta_omega; // On boundaries, phi_gamma is forced to be == phi_omega
    double x_theta = theta_omega / std::max(theta_gamma, small_real);

    double theta_velocity_omega = theta_omega * uR;
    double cR = pressureWaveCelerity(theta_omega, phi_omega);
    double lambda_in = std::max(lambdaPlus(uR, cR, x_phi, x_theta, bR), small_real);
    //== compute characteristic variables

    flux1 = 0;

    if (bc_type == "fr")
    {
        // Critical section
        flux1 = -bc_val * cR * theta_omega;
        flux2 = specificForce(theta_velocity_omega, uR, pressure_force)+
                      lambda_in * (flux1 - theta_velocity_omega);
    }
    else if (bc_type == "h")
    {
        // Prescribed h
        double h_bc = std::max(bc_val, 0.);
        double z_bc = h_bc + cell_z_tab[idx_cell][0]; // z_tab[0] identical for edges and cells at the boundary
        double theta_omega_bc = thetaWOfZ(idx_cell, z_bc);

        if (lambda_in < small_real)
        {
            flux1 = theta_velocity_omega;
            flux2 = specificForce(theta_velocity_omega, uR, pressure_force);
        }
        else
        {
            flux1 = theta_velocity_omega + lambda_in * (theta_omega_bc - theta_omega);
            flux2 = specificForce(theta_velocity_omega, uR, pressure_force) +
                          lambda_in * (flux1 - theta_velocity_omega);
        }
    }
    else if (bc_type == "z")
    {
        // Prescribed z
        double z_bc = bc_val;
        double theta_omega_bc = thetaWOfZ(idx_cell, z_bc);

        if (lambda_in < small_real)
        {
            flux1 = theta_velocity_omega;
            flux2 = specificForce(theta_velocity_omega, uR, pressure_force);
        }
        else
        {
            flux1 = theta_velocity_omega + lambda_in * (theta_omega_bc - theta_omega);
            flux2 = specificForce(theta_velocity_omega, uR, pressure_force) +
                                lambda_in * (flux1 - theta_velocity_omega);
        }
    }
    else if (bc_type == "q")
    {
        // Prescribed q
        double lambda_critic = std::exp((1./3.) * std::log(gravity * std::max(bc_val, 0.)));
        double lambda = std::max(lambda_in, lambda_critic);
        flux1 = bc_val;
        flux2 = specificForce(theta_velocity_omega, uR, pressure_force) +
                      lambda * (flux1 - theta_velocity_omega);
    }
    else
    {
        logInfo() << "BC Type " << bc_type << " not handled ";
        if (bc_type == "c")
        {
            logWarn() << "Warning 'c-type' boundary condition has been renamed as 'fr'";
        }
        if (bc_type == "f")
        {
            logWarn() << "Warning 'f-type' boundary condition has been deactivated";
        }
        logInfo() << "Edge index" << idx_edge_bnd << ":";
        logInfo() << "  x=" << edge_bnd_midpoint_x[idx_edge_bnd] << "y=" << edge_bnd_midpoint_y[idx_edge_bnd];
        logInfo() << "  Stopping...";
        return EXIT_FAILURE;
    }
    //== Routine Core
    //== get boundary condition information

    //-- compute flux of tangeantial momentum
    if (flux1 > 0)
    {
        flux3 = 0;
    }
    else
    {
        flux3 = flux1 * vR;
    }
    //== compute flux of tangeantial momentum
    return EXIT_SUCCESS;
}


bool swDdp::hypCompute(std::string const& div_corr_option, size_t n_it_div_max, std::string const& num_scheme, double& dt)
/**
 * @brief Computes the balance using fluxes to advance the solution in time
 * @param[in] div_corr_option : type of divergence correction
 * @param[in] n_it_div_max : max number of iterations allowed in the divergence correction process (if applicable)
 * @param[in] num_scheme : name of the numerical scheme
 * @param[in,out] dt Computational time step
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 *  - thetaW
 *  - thetaW_velocity_x
 *  - thetaW_velocity_y
 */
{
    double t_1 = clock();
    if (hypComputeLRStates(num_scheme) != EXIT_SUCCESS)
    {
        logError() << "Error while computing LR states in hypCompute";
        return EXIT_FAILURE;
    }
    total_t_clock_hyperb_lrstates += clock() - t_1;

    t_1 = clock();
    dt = hypTimeStep();
    total_t_clock_hyperb_timestep += clock() - t_1;
    t_1 = clock();
    if (hypComputeFluxes() != EXIT_SUCCESS)
    {
        logError() << "Error while computing fluxes in hypCompute";
        return EXIT_FAILURE;
    }
    total_t_clock_hyperb_fluxes += clock() - t_1;

    //-- Copy state variables before balance in case diverge correction is needed
    if (div_corr_option != "none")
    {
        thetaW_old = thetaW;
        thetaW_velocity_x_old = thetaW_velocity_x;
        thetaW_velocity_y_old = thetaW_velocity_y;
    }
    //== Copy state variables before balance in case diverge correction is needed

    t_1 = clock();
    if (hypComputeBalance(dt) != EXIT_SUCCESS)
        {
            logError() << "Error while computing hyperbolic balance";
            return EXIT_FAILURE;
        }
    total_t_clock_hyperb_balance += clock() - t_1;

    t_1 = clock();
    if (div_corr_option != "none")
    {
        if (hypDivCorr(div_corr_option, n_it_div_max, dt) != EXIT_SUCCESS)
        {
            logError() << "Error during divergence correction";
            return EXIT_FAILURE;
        }
    }
    total_t_clock_hyperb_divcorr += clock() - t_1;

    //---- Handle negative water depth
    if (negative_depth_detection == "active")
    {
        // Detect negative water depths and send error
        for (size_t n = 0; n < cell_nb; ++n)
        {
            if (thetaW[n] < 0)
            {
                logError() << "DDP - Negative theta" << thetaW[n]
                           << " in cell" << n;
                return EXIT_FAILURE;
            }
        }
    }
    else
    {
        // Do not detect but correctly set h=0 to avoid sqrt of negative values
        for (size_t n = 0; n < cell_nb; ++n)
        {
            if (thetaW[n] < 0.)
            {
                thetaW[n] = 0.;
                thetaW_velocity_x[n] = 0.;
                thetaW_velocity_y[n] = 0.;
            }
        }
    }
    //==== negative water depth
    return EXIT_SUCCESS;
}


/** @todo: adapt code */
bool swDdp::hypComputeBalance(const double& dt)
/**
 * @brief Computes the balance using fluxes to advance the solution in time
 * @param[in] dt Computational time step
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 *  - thetaW
 *  - thetaW_velocity_x
 *  - thetaW_velocity_y
 */
{
    //-- Prepare _old variables in case momentum dissipation is needed
    if (model_momentum_dissipation == "active")
    {
        for (size_t idx_cell = 0 ; idx_cell < cell_nb ; idx_cell++)
        {
            thetaW_old[idx_cell] = thetaW[idx_cell];
            thetaW_velocity_x_old[idx_cell] = thetaW_velocity_x[idx_cell];
            thetaW_velocity_y_old[idx_cell] = thetaW_velocity_y[idx_cell];
        }
    }
    //== Prepare _old variables in case momentum dissipation is needed
    
    //---- Boundary edges
    for (size_t idx_edge_bnd = 0; idx_edge_bnd < edge_bnd_nb; idx_edge_bnd++)
    {
        size_t idx_cell = edge_bnd_cell_idx[idx_edge_bnd];
        thetaW[idx_cell] += edge_bnd_factor_1_F1[idx_edge_bnd] * F1_bnd[idx_edge_bnd] * dt; 
        double flux2 = F2_bnd[idx_edge_bnd];
        double flux3 = F3_bnd[idx_edge_bnd];
        double source_right = s2R_bnd[idx_edge_bnd];
        thetaW_velocity_x[idx_cell] += (  edge_bnd_factor_2_F2[idx_edge_bnd] * flux2
                                      + edge_bnd_factor_2_S2[idx_edge_bnd] * source_right
                                      + edge_bnd_factor_2_F3[idx_edge_bnd] * flux3) * dt;
        thetaW_velocity_y[idx_cell] += (  edge_bnd_factor_3_F2[idx_edge_bnd] * flux2
                                      + edge_bnd_factor_3_S2[idx_edge_bnd] * source_right
                                      + edge_bnd_factor_3_F3[idx_edge_bnd] * flux3) * dt;
    }
    //==== Boundary edges
    //---- Internal edges
    for (size_t idx_edge_int = 0; idx_edge_int < edge_int_nb ; idx_edge_int++)
    {
        size_t idx_celL = edge_int_celL_idx[idx_edge_int];
        size_t idx_celR = edge_int_celR_idx[idx_edge_int];
        thetaW[idx_celL] += edge_int_factor_1_F1_L[idx_edge_int] * F1_int[idx_edge_int] * dt;
        thetaW[idx_celR] += edge_int_factor_1_F1_R[idx_edge_int] * F1_int[idx_edge_int] * dt;
        double flux2 = F2_int[idx_edge_int];
        double flux3 = F3_int[idx_edge_int];    
        double source_left = s2L_int[idx_edge_int];
        double source_right = s2R_int[idx_edge_int];
        thetaW_velocity_x[idx_celL] += (  edge_int_factor_2_F2_L[idx_edge_int] * flux2
                                      + edge_int_factor_2_S2_L[idx_edge_int] * source_left
                                      + edge_int_factor_2_F3_L[idx_edge_int] * flux3) * dt;
        thetaW_velocity_x[idx_celR] += (  edge_int_factor_2_F2_R[idx_edge_int] * flux2
                                      + edge_int_factor_2_S2_R[idx_edge_int] * source_right
                                      + edge_int_factor_2_F3_R[idx_edge_int] * flux3) * dt;
        thetaW_velocity_y[idx_celL] += (  edge_int_factor_3_F2_L[idx_edge_int] * flux2
                                      + edge_int_factor_3_S2_L[idx_edge_int] * source_left
                                      + edge_int_factor_3_F3_L[idx_edge_int] * flux3) * dt;
        thetaW_velocity_y[idx_celR] += (  edge_int_factor_3_F2_R[idx_edge_int] * flux2
                                      + edge_int_factor_3_S2_R[idx_edge_int] * source_right
                                      + edge_int_factor_3_F3_R[idx_edge_int] * flux3) * dt;
    }
    
    //==== Internal edges
    //== compute balance

    //-- Compute momentum dissipation
    if (model_momentum_dissipation == "active")
    {
        for (size_t idx_cell = 0; idx_cell < cell_nb ; idx_cell++)
        {
            if (thetaW[idx_cell] > thetaW_old[idx_cell] + small_real) 
            //-- Momentum dissipation only if h is increasing
            //-- If mu_1^2 + mu_2^2 = 0 then the flag should be set to 'inactive'
            //-- mu is a dissipation coefficient: mu=0 means that there is no dissipation
            {
                // Get the discharge increment
                double delta_q = thetaW_velocity_x[idx_cell] - thetaW_velocity_x_old[idx_cell];
                double delta_r = thetaW_velocity_y[idx_cell] - thetaW_velocity_y_old[idx_cell];
                double cos_alpha = cell_mu_cos_alpha[idx_cell];
                double sin_alpha = cell_mu_sin_alpha[idx_cell];
                double local_delta_q = cos_alpha * delta_q + sin_alpha * delta_r;
                double local_delta_r = -sin_alpha * delta_q + cos_alpha * delta_r;
                // Compute local discharge reduction
                double local_delta_q_reduction = cell_mu_1[idx_cell] * local_delta_q;
                double local_delta_r_reduction = cell_mu_2[idx_cell] * local_delta_r;
                // Compute global discharge reduction
                double delta_q_reduction = cos_alpha * local_delta_q_reduction - sin_alpha * local_delta_r_reduction;
                double delta_r_reduction = sin_alpha * local_delta_q_reduction + cos_alpha * local_delta_r_reduction;
                // Update global discharge field
                thetaW_velocity_x[idx_cell] -= delta_q_reduction;
                thetaW_velocity_y[idx_cell] -= delta_r_reduction;
            }
        }
    }
    //== Compute momentum dissipation


    //-- check balance
    /** @todo Does it really matter? */
    //---- Rounding-off error
    for (size_t n = 0; n < cell_nb; ++n)
    {
        //if (std::fabs(thetaW[n]) < small_real) thetaW[n] = 0.;
        //if (std::fabs(thetaW_velocity_x[n]) < small_real) thetaW_velocity_x[n] = 0.;
        //if (std::fabs(thetaW_velocity_y[n]) < small_real) thetaW_velocity_y[n] = 0.;
    }
    //==== Rounding-off error
    //== check balance
    return EXIT_SUCCESS;
}


bool swDdp::hypComputeFluxes(void)
/**
 * @brief Computes the fluxes by solving the Riemann problems
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 *  - all fluxes
 */
{
    //-- Boundary edges
    for (swBoundaryConditionTimeSeries* bc :bc_timeseries_list)
    {
        bc->computeValBC(t);
    }
    for (size_t idx_edge_bnd = 0; idx_edge_bnd < edge_bnd_nb; idx_edge_bnd++)
    {
        double zR = z_riemann_bnd[idx_edge_bnd];
        double uR_rotated = 0.0;
        double vR_rotated = 0.0;
        if (theta_riemann_bnd[idx_edge_bnd] > small_real)
        {
            uR_rotated = theta_u_riemann_bnd[idx_edge_bnd] / theta_riemann_bnd[idx_edge_bnd];//q_riemann_bnd[idx_edge_bnd];
            vR_rotated = theta_v_riemann_bnd[idx_edge_bnd] / theta_riemann_bnd[idx_edge_bnd];//r_riemann_bnd[idx_edge_bnd];
        }        
        double flux1 = 0.0;
        double flux2 = 0.0;
        double flux3 = 0.0;
        double source_right = 0.0;

        // Get BC type and value
        double bc_val = edge_bnd_bc_pointer[idx_edge_bnd]->getBCVal();
        std::string bc_type = edge_bnd_bc_pointer[idx_edge_bnd]->getBCType();
        size_t idx_cell = edge_bnd_cell_idx[idx_edge_bnd];
        if (hypBoundaryFlux(idx_edge_bnd, boussinesq[idx_cell], cell_z_b[idx_cell],
                         zR, uR_rotated, vR_rotated, 
                         lambda_plus_bnd[idx_edge_bnd], bc_type, bc_val,
                         flux1, flux2, flux3, source_right) != EXIT_SUCCESS)
        {
            logError() << "Error while computing boundary fluxes";
            return EXIT_FAILURE;
        }
    
        F1_bnd[idx_edge_bnd] = flux1;
        F2_bnd[idx_edge_bnd] = flux2;
        F3_bnd[idx_edge_bnd] = flux3;
        s2R_bnd[idx_edge_bnd] = source_right;        
    }
    //=== Boundary edges

    //-- Internal edges
    for (size_t idx_edge_int = 0 ; idx_edge_int < edge_int_nb ; idx_edge_int++)
    {
        size_t celL_idx = edge_int_celL_idx[idx_edge_int];
        size_t celR_idx = edge_int_celR_idx[idx_edge_int];
        double bL = boussinesq[celL_idx];
        double bR = boussinesq[celR_idx];
        double zbL = cell_z_b[edge_int_celL_idx[idx_edge_int]];
        double zbR = cell_z_b[edge_int_celR_idx[idx_edge_int]];
        double zL = z_riemann_int[idx_edge_int][0];
        double zR = z_riemann_int[idx_edge_int][1];
        double hL = zL - zbL;
        double hR = zR - zbR;
        double uL = 0.0;
        double vL = 0.0;
        if (theta_riemann_int[idx_edge_int][0]>small_real)
        {
            uL = theta_u_riemann_int[idx_edge_int][0] / theta_riemann_int[idx_edge_int][0];
            vL = theta_v_riemann_int[idx_edge_int][0] / theta_riemann_int[idx_edge_int][0];
        }
        double uR = 0.0;
        double vR = 0.0;
        if (theta_riemann_int[idx_edge_int][1]>small_real)
        {
            uR = theta_u_riemann_int[idx_edge_int][1] / theta_riemann_int[idx_edge_int][1];
            vR = theta_v_riemann_int[idx_edge_int][1] / theta_riemann_int[idx_edge_int][1];
        }
        // Compute fluxes only if sufficient depth

        if ((hL > h_min) || (hR > h_min))
        {
            double flux1 = 0;
            double flux2 = 0;
            double flux3 = 0;
            double source_left = 0;
            double source_right = 0;

            if (topo_discretization == "discontinuous")
            {
                if (hypHllc(idx_edge_int, bL, bR, zL, uL, vL, zR, uR, vR,
                        flux1, flux2, flux3,
                        source_left, source_right) != EXIT_SUCCESS)
                {
                    logError() << "Error in HLLC";
                    return EXIT_FAILURE;
                }
            }
            /** @todo Blending with DDP has been removed for now...
            else if (topo_discretization == "blending")
            {
                hypHllcBlending(bL, bR, phiG, phiL, phiR, zbL, zbR,
                                zL, qL, rL, zR, qR, rR,
                                lambda_minus_LR[0], lambda_minus_LR[1],
                                lambda_plus_LR[0], lambda_plus_LR[1],
                                flux1, flux2, flux3,
                                source_left, source_right);
            }
            */ 
            else
            {
                logError() << "'topo_discretization' parameter value not recognized" << topo_discretization;
                return EXIT_FAILURE;
            }

            if (std::isnan(flux1) == true)
            {
                logError() << "DDP - Mass flux is Not a number for internal interface " << idx_edge_int;
                logInfo() << "DEBUG INTERFACE" << edge_int_idx[idx_edge_int]
                    << "nx=" << edge_int_norm_x[idx_edge_int]
                    << "ny=" << edge_int_norm_y[idx_edge_int];
                logInfo() << "  Left-hand cell" << edge_int_celL_idx[idx_edge_int]
                    << "zb=" << zbL
                    << "z=" << zL ;
                logInfo() << "  Right-hand cell " << edge_int_celR_idx[idx_edge_int]
                    << "zb=" << zbR
                    << "z=" << zR ;
                logInfo() << "  flux1=" << flux1
                    << "flux2=" << flux2
                    << "flux3=" << flux3
                    << "S2L=" << source_left
                    << "S2R=" << source_right;
                return EXIT_FAILURE;
            }
            F1_int[idx_edge_int] = flux1;
            F2_int[idx_edge_int] = flux2;
            F3_int[idx_edge_int] = flux3;
            s2L_int[idx_edge_int] = source_left;
            s2R_int[idx_edge_int] = source_right;
        }
    }
    //== Internal edges
    return EXIT_SUCCESS;
}


bool swDdp::hypComputeLRStates(std::string const& num_scheme)
/**
 * @brief Defines the local Riemann problems for the given numerical scheme
 * @param[in] num_scheme Numerical scheme used (e.g. godunov, muscl-evr)
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 *  @todo to be checked
 */
{
    if (num_scheme == "godunov")
    {
        if (hypGodunovComputeLRStates() != EXIT_SUCCESS)
        {
            logError() << "Error while computing Godunaov LR states";
            return EXIT_FAILURE;
        }
    }
    else if (num_scheme == "muscl_evr")
    {
        //hypMusclEvrComputeLRStates();
        logError() << "MUSCL EVR numerical scheme is not implemented : " << num_scheme;
        return EXIT_FAILURE;
    }
    else
    {
        logError() << "Numerical scheme is not implemented : " << num_scheme;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


double swDdp::hypTimeStep(void)
/**
 * @brief Computes the CFL-permissible computational timestep for the DDP model
 * @return dt : computational time step
 */
{
    double dt = dt_max;
    std::vector<double> swept(cell_nb, 0.);

    // --- Boundary edges
    for (size_t idx_edge_bnd = 0; idx_edge_bnd < edge_bnd_nb; idx_edge_bnd++)
    {
        size_t id_cell = edge_bnd_cell_idx[idx_edge_bnd];
        double lambda = 0.;
        double h1 = z_riemann_bnd[idx_edge_bnd] - cell_z_tab[id_cell][0];
        if (h1 > h_min)
        {
            lambda = lambda_max_bnd[idx_edge_bnd];
            swept[id_cell] += lambda * edge_bnd_length[idx_edge_bnd];
        }
    }
    // === Boundary edges

    // --- Internal edges
    for (size_t idx_edge_int = 0; idx_edge_int < edge_int_nb; idx_edge_int++)
    {
        double lambda = 0.;
        
        // --- Left-hand side
        size_t idx_celL = edge_int_celL_idx[idx_edge_int];
        double h1 = z_riemann_int[idx_edge_int][0] - cell_z_tab[idx_celL][0];
        if (h1 > h_min)
        {
            lambda = lambda_max_int[idx_edge_int];
        }
        swept[idx_celL] += lambda * edge_int_length[idx_edge_int];
        // === Left-hand side

        // --- Right-hand side
        size_t idx_celR = edge_int_celR_idx[idx_edge_int];
        h1 = z_riemann_int[idx_edge_int][1] - cell_z_tab[idx_celR][0];
        if (h1 > h_min)
        {
            lambda = lambda_max_int[idx_edge_int];
        }
        // === Right-hand side
        swept[idx_celR] += lambda * edge_int_length[idx_edge_int];
    }
    // === Internal edges

    // Cells : determine the maximum Dt
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        double xxx = cfl_max * (cell_area[cell_id[idx_cell]] /
            std::max(swept[cell_id[idx_cell]], small_real));
        if (xxx < dt)
        {
            dt = xxx;
            if (dt < 0.)
            {
                logError() << "Negative Time step in CFL limiting!  deltat :" << dt << "in cell:" << cell_id[idx_cell];
            }
        }
    }

    //-- Map storage-related time step limitation
    store_map = false;
    if (dt_map <= -small_real) 
    {
        // Negative dt means maps will be stored every time step
        store_map = true;
    }
    else if (dt_map > small_real)
    {
        int t_now_over_dt_int = (int)t / dt_map;
        int t_later_over_dt_int = (int)(t + dt) / dt_map;
        if (t_later_over_dt_int > t_now_over_dt_int) /* Storing is mandatory */
        {
            dt = t_later_over_dt_int * dt_map - t; /* Revise dt so as to reach exactly the nearest storage time */
            store_map = true;
        }
    }
    if ((t + dt) > t_max)
    {
        dt = t_max - t;
    }
    //== Map storage-related time step limitation

    //-- Determine whether time series should be stored
    store_time_series = false;
    if (dt_time_series <= -small_real)
    {
        // Negative dt_time_series time series will be stored every time step
        store_time_series = true;
    }
    else if (dt_time_series > small_real)
    {
        int t_now_over_dt_int = (int)t / dt_time_series;
        int t_later_over_dt_int = (int)(t + dt) / dt_time_series;
        if (t_later_over_dt_int > t_now_over_dt_int) /* Storing is mandatory */
        {
            store_time_series = true;
        }
    }
    //== Determine whether time series should be stored

    return dt;
}


/** @todo: adapt code */
bool swDdp::hypDivCorr(std::string const& div_corr_option, size_t n_it_div_max, double& dt)
/**
 * @brief Corrects fluxes or time step from divergence considerations
 * @param[in] div_corr_option The desired option for divergence correction
 * @param[in,out] dt Computational time step
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 */
{
    if (div_corr_option == "limit_timestep")
    {
        if (hypDivCorrTimeStep() != EXIT_SUCCESS)
        {
            logError() << "Error during time-step modification for divergence correction";
            return EXIT_FAILURE;
        }
    }
    else if (div_corr_option == "limit_flux")
    {
        if (hypDivCorrFluxes(n_it_div_max) != EXIT_SUCCESS)
        {
            logError() << "Error during fluxes modification for divergence correction";
            return EXIT_FAILURE;
        }
    }
    else
    {
        logError() << "Divergence correction option is not handled: " << div_corr_option;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swDdp::hypDivCorrFluxes(size_t n_it_div_max)
/**
 * @brief compute the divergence correction (with the fluxes limitation) to avoid overdrying
 * @param[in] n_it_div_max : max number of iterations in the process
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes :
 *  - all fluxes
 */
{
    bool correct = false;
    size_t nbiter = 0;

    while (correct == false)
    {
        size_t ndry = 0;
        for (size_t idx_cell = 0 ; idx_cell < cell_nb ; idx_cell++)
        {
            //-- For each cell, compute the flux factor that avoids drying and revise fluxes accordingly
            double flux_ratio = 1.0;
            if (thetaW[idx_cell] < 0.0)
            {
                ndry++;
                flux_ratio = thetaW_old[idx_cell] / (thetaW_old[idx_cell] - thetaW[idx_cell]);
                std::vector<size_t> id_edge_bnd = cell_edges_bnd_idx[idx_cell];
                for (size_t idx_edge_bnd = 0; idx_edge_bnd < id_edge_bnd.size(); idx_edge_bnd++)
                {
                    F1_bnd[idx_edge_bnd] *= flux_ratio;
                    F2_bnd[idx_edge_bnd] *= flux_ratio;
                    F3_bnd[idx_edge_bnd] *= flux_ratio;
                    s2R_bnd[idx_edge_bnd] *= flux_ratio;
                }
                std::vector<size_t> id_edge_int = cell_edges_int_idx[idx_cell];
                for (size_t idx_edge_int = 0; idx_edge_int < id_edge_int.size(); idx_edge_int++)
                {
                    F1_int[idx_edge_int] *= flux_ratio;
                    F2_int[idx_edge_int] *= flux_ratio;
                    F3_int[idx_edge_int] *= flux_ratio;
                    s2L_int[idx_edge_int] *= flux_ratio;
                    s2R_int[idx_edge_int] *= flux_ratio;
                }
            }
            //== For each cell, compute the flux factor that avoids drying and revise fluxes accordingly
        }
        ++nbiter;
        //-- Recompute balance
        for (size_t idx_edge_int = 0; idx_edge_int < edge_int_nb; idx_edge_int++)
        {
            size_t idx_celL = edge_int_celL_idx[idx_edge_int];
            size_t idx_celR = edge_int_celR_idx[idx_edge_int];
            thetaW[idx_celL] = thetaW_old[idx_celL] + edge_int_factor_1_F1_L[idx_edge_int] * F1_int[idx_edge_int] * dt;
            thetaW[idx_celR] = thetaW_old[idx_celR] + edge_int_factor_1_F1_R[idx_edge_int] * F1_int[idx_edge_int] * dt;
            thetaW_velocity_x[idx_celL] = thetaW_velocity_x_old[idx_celL] + (  edge_int_factor_2_F2_L[idx_edge_int] * F2_int[idx_edge_int]
                                        + edge_int_factor_2_S2_L[idx_edge_int] * s2L_int[idx_edge_int]
                                        + edge_int_factor_2_F3_L[idx_edge_int] * F3_int[idx_edge_int]) * dt;
            thetaW_velocity_x[idx_celR] = thetaW_velocity_x_old[idx_celR] + (  edge_int_factor_2_F2_R[idx_edge_int] * F2_int[idx_edge_int]
                                        + edge_int_factor_2_S2_R[idx_edge_int] * s2R_int[idx_edge_int]
                                        + edge_int_factor_2_F3_R[idx_edge_int] * F3_int[idx_edge_int]) * dt;
            thetaW_velocity_y[idx_celL] = thetaW_velocity_y_old[idx_celL] + (  edge_int_factor_3_F2_L[idx_edge_int] * F2_int[idx_edge_int]
                                        + edge_int_factor_3_S2_L[idx_edge_int] * s2L_int[idx_edge_int]
                                        + edge_int_factor_3_F3_L[idx_edge_int] * F3_int[idx_edge_int]) * dt;
            thetaW_velocity_y[idx_celR] = thetaW_velocity_y_old[idx_celR] +  (  edge_int_factor_3_F2_R[idx_edge_int] * F2_int[idx_edge_int]
                                        + edge_int_factor_3_S2_R[idx_edge_int] * s2R_int[idx_edge_int]
                                        + edge_int_factor_3_F3_R[idx_edge_int] * F3_int[idx_edge_int]) * dt;
        }
        //== Recompute balance
        if ((ndry == 0) || (nbiter > n_it_div_max))
        {
            correct = true;
        }
    }
    return EXIT_SUCCESS;
}


bool swDdp::hypDivCorrTimeStep(void)
/**
 * @brief Corrects time step to prevent overdrying
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes :
 *  - dt
 */
{
    double dt_initial = dt;
    double dt_avoid_drying = dt;
    bool some_cells_will_be_dry = false;
    //-- Check whether some cells will become dry, correct dt as appropriate
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        size_t id_cell = cell_id[idx_cell];
        if (thetaW[id_cell] <= 0.0)
        {
            double dt_avoid_drying_cell = dt * (thetaW_old[id_cell] - thetaW[id_cell]) / thetaW[id_cell];
            dt_avoid_drying = std::min(dt_avoid_drying_cell, dt_avoid_drying);
            some_cells_will_be_dry = true;
        }
    }
    dt = std::min(dt, dt_avoid_drying);
    //== Check whether some cells will become dry, correct dt as appropriate

    //-- Recompute variables at next time level if needed
    if (some_cells_will_be_dry)
    {
        double dt_ratio = dt / dt_initial;
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            size_t id_cell = cell_id[idx_cell];
            thetaW[id_cell] = thetaW_old[id_cell] * (1.0 - dt_ratio) + thetaW[id_cell] * dt_ratio;
            thetaW_velocity_x[id_cell] = thetaW_velocity_x_old[id_cell] * (1.0 - dt_ratio) + thetaW_velocity_x[id_cell] * dt_ratio;
            thetaW_velocity_y[id_cell] = thetaW_velocity_y_old[id_cell] * (1.0 - dt_ratio) + thetaW_velocity_y[id_cell] * dt_ratio;
        }
        //== Recompute variables at next time level if needed
    }
    return EXIT_SUCCESS;
}


/** @todo: adapt code */
double swDdp::hypHCriticalPhiq(double const& q, double const& phiG, double const& PhiW)
/**
 * @brief swDdpHyp::hCritical
 * @param q: unit discharge
 * @param phiG: interface porosity
 * @param PhiW:  cell porosity
 * @return
 */
{
    double h_critical = q * PhiW / std::max(phiG, small_real);
    //-- @todo: remove once proven ok:  h_critical = h_critical * h_critical / m_param->g;
    h_critical = h_critical * h_critical / gravity;
    h_critical = std::exp((1./3.) * std::log(h_critical));
    return h_critical;
}


/** @todo: adapt code */
double swDdp::hypHCriticalPhiWq(double const& phiwq, double const& phiG)
/**
 * @brief swDdpHyp::hCritical
 * @param q: unit discharge
 * @param phiG: interface porosity
 * @param PhiW:  cell porosity
 * @return
 */{
    double h_critical = phiwq / std::max(phiG, small_real);
    //-- @todo: remove once proven ok:  h_critical = h_critical * h_critical / m_param->g;
    h_critical = h_critical * h_critical / gravity;
    h_critical = std::exp((1./3.) * std::log(h_critical));
    return h_critical;
}


/** @todo Shall we pass lamndas as parameters? */
bool swDdp::hypHllc(const size_t& idx_edge_int, double const& bL, double const& bR,
                         double const& zL, double const& uL, double const& vL,
                         double const& zR, double const& uR, double const& vR, 
                         double& flux1, double& flux2, double& flux3,
                         double& source_left, double& source_right)
/**
 * @brief computes the fluxes using the classical HLLC formulation in Guinot & Soares-Frazao (2006)
 * @param[in] bL : Boussinesq coefficient in the cell on the LHS of the interface
 * @param[in] bR : Boussinesq coefficient in the cell on the RHS of the interface
 * @param[in] zL  : free surface elevation in the cell on the LHS of the interface
 * @param[in] uL  : normal velocity in the cell on the LHS of the interface
 * @param[in] vL  : tangent uvelocity in the cell on the LHS of the interface
 * @param[in] zR  : free surface elevation in the cell on the RHS of the interface
 * @param[in] uR  : normal velocity in the cell on the RHS of the interface
 * @param[in] vR  : tangent velocity in the cell on the RHS of the interface
 * @param[out] flux1 : mass flux in the direction normal to the interface
 * @param[out] flux2 : momentum flux in the direction normal to the interface (specific force per unit width)
 * @param[out] flux3 : momentum flux in the direction tangent to the interface
 * @param[out] source_left: geometric source term on the LHS of the interface
 * @param[out] source_right: geometric source term on the RHS of the interface
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 */
{
    size_t idx_celL = edge_int_celL_idx[idx_edge_int];
    size_t idx_celR = edge_int_celR_idx[idx_edge_int];

    double phi_omega_L = phiWOfZ(idx_celL, zL);
    double phi_omega_R = phiWOfZ(idx_celR, zR);
    double theta_omega_L = thetaWOfZ(idx_celL, zL);
    double theta_omega_R = thetaWOfZ(idx_celR, zR);
    //== initialize parameters
    double theta_velocity_omega_L = uL * theta_omega_L;
    double theta_velocity_omega_R = uR * theta_omega_R;


    //-- compute flow variables - Gamma
    double theta_gamma_L = thetaGIntOfZ(idx_edge_int, zL);
    double theta_gamma_R = thetaGIntOfZ(idx_edge_int, zR);
    double theta_gamma_hmin = thetaG_int_hmin[idx_edge_int];
    double theta_velocity_gamma_L = theta_velocity_omega_L * std::min(theta_gamma_L / theta_gamma_hmin, 1.0);
    double theta_velocity_gamma_R = theta_velocity_omega_R * std::min(theta_gamma_R / theta_gamma_hmin, 1.0);

    double velocity_gamma_L = theta_velocity_gamma_L/std::max(theta_gamma_L,small_real);
    double velocity_gamma_R = theta_velocity_gamma_R/std::max(theta_gamma_R,small_real);

    double z_gamma_min = edge_int_z_tab[idx_edge_int][0];
    if (zL - z_gamma_min < h_min)
    {
        double umax = std::sqrt(uL * uL + 2. * gravity * std::max(zL-z_gamma_min,small_real));
        velocity_gamma_L = std::min(umax,std::max(velocity_gamma_L,-umax));
        theta_velocity_gamma_L = velocity_gamma_L * theta_gamma_L;
    }
    if (zR - z_gamma_min < h_min)
    {
        double umax = std::sqrt(uR * uR + 2. * gravity * std::max(zR-z_gamma_min,small_real));
        velocity_gamma_R = std::min(umax,std::max(velocity_gamma_R,-umax));
        theta_velocity_gamma_R = velocity_gamma_R * theta_gamma_R;
    }
    //== compute flow variables - Gamma
    
    //-- Pressure forces (water-water and solid-water)
    double Fp_L = 0;
    double Fp_R = 0;
    double Fb_L = 0;
    double Fb_R = 0;
    if (computePressureForceOmegaOntoGammaOfZ(idx_edge_int, zL, zR, Fp_L, Fp_R, Fb_L, Fb_R) != EXIT_SUCCESS)
    {
        logError() << "Error while computing pressure force from Omega to Gamma";
        return EXIT_FAILURE;
    }
    //== Pressure forces (water-water and solid-water)

    //-- Compute fluxes
    //---- Left and right states for the Riemann problem
    double M_L = specificForce(theta_velocity_gamma_L, velocity_gamma_L, Fp_L);
    double M_R = specificForce(theta_velocity_gamma_R, velocity_gamma_R, Fp_R);
    double zzL = std::max(zL,z_gamma_min);
    double zzR = std::max(zR,z_gamma_min);
    //==== Left and right states for the Riemann problem

    double theta_RL = thetaWOfZ(idx_celR, zL); 
    double theta_LR = thetaWOfZ(idx_celL, zR);
    // soit L peut passer au dessus de l'interface, soit R => max entre les deux, pas min
    // mais il faut qu'il y ait assez de place dans la cellule d'à côté => min entre interface et cellule
    //double min_theta = std::min(theta_RL, std::min(theta_LR, std::max(theta_gamma_L, theta_gamma_R)));
    if ((theta_gamma_L > theta_RL) ||
        (theta_gamma_R > theta_LR) ||
        (theta_gamma_L > theta_omega_L) ||
        (theta_gamma_R > theta_omega_R))
    {
        warning_thetaG_gt_thetaW += 1;
    }

    double min_theta = std::max(std::min(theta_RL,theta_gamma_L), std::min(theta_LR,theta_gamma_R));
    //Todo: verifier que min(theta_RL,theta_gamma_L) est fait au moment de la construction des phi gamma
    // si c'est ok: min_theta = std::max(theta_gamma_L, theta_gamma_R);
    flux1 = 0;
    if (((theta_omega_L > small_real) && (theta_omega_R > small_real)) ||
         (min_theta >= small_real))
    {
        //------ There is water in both cells or one of the cells can fill the other
        double cL = pressureWaveCelerity(theta_omega_L, phi_omega_L);
        double cR = pressureWaveCelerity(theta_omega_R, phi_omega_R);

        double phi_gamma_L = phiGIntOfZ(idx_edge_int, zL);
        double phi_gamma_R = phiGIntOfZ(idx_edge_int, zR);
        double chi_theta_L = theta_omega_L / std::max(theta_gamma_L, theta_gamma_hmin);
        double chi_theta_R = theta_omega_R / std::max(theta_gamma_R, theta_gamma_hmin);
        double chi_phi_L = phi_omega_L / std::max(phi_gamma_L,small_real);
        double chi_phi_R = phi_omega_R / std::max(phi_gamma_R,small_real);
        chi_phi_L = std::min(chi_phi_max, chi_phi_L);
        chi_phi_R = std::min(chi_phi_max, chi_phi_R);

        double L1L = lambdaMinus(uL, cL, chi_phi_L, chi_theta_L, bL);
        double L2L = lambdaPlus(uL, cL, chi_phi_L, chi_theta_L, bL);
        double L1R = lambdaMinus(uR, cR, chi_phi_R, chi_theta_R, bR);
        double L2R = lambdaPlus(uR, cR, chi_phi_R, chi_theta_R, bR);
 
        double L1 = std::min({L1L, L1R, 0.});   // return the min of the list {}
        double L2 = std::max({L2L, L2R, 0.});

        flux1 = (L2 * phi_omega_R * theta_velocity_gamma_L -
                             L1 * phi_omega_L * theta_velocity_gamma_R +
                             L1 * L2 * std::min(phi_omega_L, phi_omega_R) * std::min(phi_gamma_L, phi_gamma_R) * (zzR - zzL)) /
                             std::max(phi_omega_R * L2 - phi_omega_L * L1, small_real);
        flux2 = (L2 * M_L - L1 * M_R +
                             L1 * L2 * (theta_omega_R * uR - theta_omega_L * uL)) /
                            std::max(L2 - L1, small_real);
        if (flux1 < 0)
        {
            //-------- Mass fluxes from right to left
            flux3 = flux1 * vR;
            //======== Mass fluxes from right to left
        }
        else
        {
            //-------- Mass fluxes from left to right
            flux3 = flux1 * vL;
            //======== Mass fluxes from left to right
        }
        source_left = -Fb_L;
        source_right = Fb_R;
        if ((theta_gamma_L < small_real) &&
            (theta_gamma_R < small_real))
        {
            // water in cell but impervious interface
            flux1 = 0.;
            flux2 = 0.;
            flux3 = 0.;
        }
        //====== There is water in both cells or one of the cells can fill the other
    }
    else if ((theta_omega_L < small_real) &&
             (min_theta < small_real))
    {
        /* The cell L is empty and the water level in the cell R is below zbL
         * => filling is impossible*/
        flux1 = 0.;
        flux2 = 0.;
        flux3 = 0.;
        source_left = 0.;
        source_right = Fb_R;

    }
    else if ((theta_omega_R < small_real) &&
             (min_theta < small_real))
    {
        /* The cell R is empty and the water level in the cell L is below zbR
         * => filling is impossible*/
        flux1 = 0.;
        flux2 = 0.;
        flux3 = 0.;
        source_left = - Fb_L;
        source_right = 0.;
    }
    else
    {
        logError() << "Impossible case";
        logError() << small_real;
        logError() << "id_edge_int = " << idx_edge_int;
        logError() << "Left cell = " << idx_celL;
        logError() << "Left Z_TAB = " << cell_z_tab[idx_celL];
        logError() << "Right cell = " << idx_celR;
        logError() << "Right Z_TAB = " << cell_z_tab[idx_celR];
        logError() << "zL=" << zL << "\ttheta_omega_L=" << theta_omega_L << "\ttheta_gamma_L=" << theta_gamma_L << "\ttheta_RL" << theta_RL;
        logError() << "zR=" << zR << "\ttheta_omega_R=" << theta_omega_R << "\ttheta_gamma_R=" << theta_gamma_R << "\ttheta_LR" << theta_LR;
        logError() << "min_theta=" << min_theta;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}