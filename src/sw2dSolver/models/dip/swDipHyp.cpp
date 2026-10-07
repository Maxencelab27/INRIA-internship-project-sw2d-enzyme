#include <algorithm>
#include <cmath>
#include <vector>

#include "./swDip.h"
#include "../common/swSource.h"
#include "../common/swBoundaryConditionTimeSeries.h"
#include "../../../common/swBoundaryCondition.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"
#include "../../../common/param.h"

void swDip::hypBoundaryFlux(const double& bR, const double& phiG,const double& phiR, const double& zbR, const double& zR, const double& qR, const double& rR,
                           const double& lambda_plus_R, 
                           std::string const& bc_type, const double& bc_val,
                           double& flux1, double& flux2, double& flux3, double& source_right)
/**
 * @brief Computes fluxes at boundary edges for the DIP models
 * @param[in] bR  : Boussinesq coefficient on the RHS of the interface
 * @param[in] phiG  : connectivity porosity at the boundary edge
 * @param[in] phiR  : storage porosity in the cell
 * @param[in] zbR : bottom level in the cell on the RHS of the interface
 * @param[in] zR  : free surface elevation in the cell on the RHS of the interface
 * @param[in] qR  : normal unit discharge in the cell on the RHS of the interface
 * @param[in] rR  : tangent unit discharge in the cell on the RHS of the interface
 * @param[in] lambda_plus_R : lambda+ value in the cell on the RHS of the interface
 * @param[in] bc_type : type of boundary condition on the boundary interface ('h', 'q', 'fr', etc.)
 * @param[in] bc_val : value for the boundary condition of type bc_type
 * @param[out] flux1 : mass flux in the direction normal to the interface
 * @param[out] flux2 : momentum flux in the direction normal to the interface (specific force per unit width)
 * @param[out] flux3 : momentum flux in the direction tangent to the interface
 * @param[out] source_right: geometric source term on the RHS of the interface
 */
    /* Formulae are designed for a domain on the right-hand side of the boundary */
{
    //-- recompute hydrodynamic values in the computational cell
    double hR = zR - zbR;
    double q1 = qR;
    double r1 = rR;

    if (hR < h_min)
    {
        q1 = 0.;
        r1 = 0.;
        hR = 0.;
    }
    //== recompute hydrodynamic values in the computational cell

    //-- compute characteristic variables
    double cR = std::sqrt(gravity * std::max(hR, small_real));
    double uR = q1 / std::max(hR, small_real);
    double L1 = std::max(0.0, lambda_plus_R);
    double flux2R = specificForce(bR, phiG, phiR, hR, uR);
    //== compute characteristic variables

    //-- Routine Core
    if (bc_type == "fr")
    {
        // Critical section
        double c = (uR - 2. * cR) / (- bc_val - 2.);
        c = std::max(c, 0.);
        double u = - bc_val * c;
        double h = c * c / gravity;
        flux1 = phiR * u * h;
        flux2 = flux2R + (flux1 - q1 * phiR) * L1;
    }
    else if (bc_type == "h")
    {
        // Prescribed h
        double h = std::max(bc_val, 0.);
        flux1 = phiR * (q1 + L1 * (h - hR));
        flux2 = flux2R + (flux1 - q1 * phiR) * L1;
        /*if (debug_bc == 12)
        {
            logInfo() << "     BC2 - TypBC  ValBC " << typ << bc_val;
            logInfo() << "           zb           " << zbR;
            logInfo() << "           h1, q1       " << hR << q1;
            logInfo() << "           L1, h        " << L1 << h;
            logInfo() << "           flux1, flux2       " << flux1 << flux2;
        }*/
    }
    else if (bc_type == "z")
    {
        // Prescribed z
        double h = std::max(bc_val - zbR, 0.);
        flux1 = phiR * (q1 + L1 * (h - hR));
        flux2 = flux2R + (flux1 - q1 * phiR) * L1;
        /*if (debug_bc == 12)
        {
            logInfo() << "     BC2 - TypBC  ValBC " << typ << bc_val;
            logInfo() << "           zb           " << zbR;
            logInfo() << "           h1, q1       " << hR << q1;
            logInfo() << "           L1, h        " << L1 << h;
            logInfo() << "           flux1, flux2       " << flux1 << flux2;
        }*/
    }
    else if (bc_type == "phiWq")
    {
        // Prescribed q
        double h_critical = hypHCriticalPhiWq(bc_val, phiG);
        if (h_critical >= hR)
        {
            // water depth in the domain is low => discharge entering with the critical depth
            flux1 = phiG * 1.5 * std::sqrt(gravity * h_critical) * h_critical;
            flux2 = flux2R +  (flux1 - phiR * q1) * L1;
        }
        else
        {
            flux1 = bc_val;
            flux2 = flux2R +  (flux1 - phiR * q1) * L1;
        }
    }
    else if (bc_type == "q")
    {
        // Prescribed q
        double h_critical = hypHCriticalPhiq(bc_val, phiG, phiR);
        if (h_critical >= hR)
        {
            // water depth in the domain is low => discharge entering with the critical depth
            flux1 = phiG * 1.5 * std::sqrt(gravity * h_critical) * h_critical;
            flux2 = flux2R +  (flux1 - phiR * q1) * L1;
        }
        else
        {
            flux1 = phiR * bc_val;
            flux2 = flux2R +  (flux1 - phiR * q1) * L1;
        }
    }
    else
    {
        logError() << "BC Type " << bc_type << " not handled ";
        logFatal(__FILE__, __func__, __LINE__);
        
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
        double vR = r1 / std::max(hR, small_real);
        flux3 = flux1 * vR;
    }
    //== compute flux of tangeantial momentum

    //-- Add reaction of solid boundary
    flux2 += std::max(phiR - phiG, 0.) * 0.5 * gravity * hR * hR;
    /** @todo Check validity of source code line below, this was apparently forgotten in the old code */
    source_right = 0.0;
    //== Add reaction of solid boundary
}


void swDip::hypCompute(std::string const& div_corr_option, size_t n_it_div_max, std::string const& num_scheme, double& dt)
/**
 * @brief Computes the balance using fluxes to advance the solution in time
 * @param[in] div_corr_option : type of divergence correction
 * @param[in] n_it_div_max : max number of iterations allowed in the divergence correction process (if applicable)
 * @param[in] num_scheme : name of the numerical scheme
 * @param[in,out] dt Computational time step
 * @return nothing
 * @par Modified attributes:
 *  - water_height
 *  - unit_discharge_x
 *  - unit_discharge_y
 */
{
    double t_1 = clock();
    hypComputeLRStates(num_scheme);
    total_t_clock_hyperb_lrstates += clock() - t_1;

    t_1 = clock();
    dt = hypComputeTimeStep();
    total_t_clock_hyperb_timestep += clock() - t_1;

    t_1 = clock();
    hypComputeFluxes();
    total_t_clock_hyperb_fluxes += clock() - t_1;

    //-- Copy state variables before balance in case diverge correction is needed
    if (div_corr_option != "none")
    {
        water_depth_old = water_depth;
        unit_discharge_x_old = unit_discharge_x;
        unit_discharge_y_old = unit_discharge_y;
    }
    //== Copy state variables before balance in case diverge correction is needed

    t_1 = clock();
    hypComputeBalance(dt);
    total_t_clock_hyperb_balance += clock() - t_1;

    t_1 = clock();
    if (div_corr_option != "none")
    {
        hypDivCorr(div_corr_option, n_it_div_max, dt);
    }
    total_t_clock_hyperb_divcorr += clock() - t_1;

    //---- Handle negative water depth
    if (negative_depth_detection == "active")
    {
        // Detect negetive water depths and send error
        for (size_t n = 0; n < cell_nb; ++n)
        {
            if (water_depth[n] < 0)
            {
                logError() << "DIP - Negative water depth" << water_depth[n]
                           << " in cell" << n;
                logFatal(__FILE__, __func__, __LINE__);
            }
        }
    }
    else
    {
        // Do not detect but correctly set h=0 to avoid sqrt of negative values
        for (size_t n = 0; n < cell_nb; ++n)
        {
            if (water_depth[n] < 0.)
            {
                water_depth[n] = 0.;
                unit_discharge_x[n] = 0.;
                unit_discharge_y[n] = 0.;
            }
        }
    }
    //==== negative water depth
}


void swDip::hypComputeBalance(const double& dt)
/**
 * @brief Computes the balance using fluxes to advance the solution in time
 * @param[in] dt Computational time step
 * @return nothing
 * @par Modified attributes:
 *  - water_height
 *  - unit_discharge_x
 *  - unit_discharge_y
 */
{
    //-- Prepare _old variables in case momentum dissipation is needed
    if (model_momentum_dissipation == "active")
    {
        for (size_t idx_cell = 0 ; idx_cell < cell_nb ; idx_cell++)
        {
            water_depth_old[idx_cell] = water_depth[idx_cell];
            unit_discharge_x_old[idx_cell] = unit_discharge_x[idx_cell];
            unit_discharge_y_old[idx_cell] = unit_discharge_y[idx_cell];
        }
    }
    //== Prepare _old variables in case momentum dissipation is needed

    //---- Boundary edges
    for (size_t edge_idx = 0; edge_idx < edge_bnd_nb; edge_idx++)
    {
        size_t id_cell = edge_bnd_cell_idx[edge_idx];
        water_depth[id_cell] += edge_bnd_factor_1_F1[edge_idx] * F1_bnd[edge_idx] * dt;
        unit_discharge_x[id_cell] += (  edge_bnd_factor_2_F2[edge_idx] * F2_bnd[edge_idx]
                                      + edge_bnd_factor_2_S2[edge_idx] * s2R_bnd[edge_idx]
                                      + edge_bnd_factor_2_F3[edge_idx] * F3_bnd[edge_idx]) * dt;
        unit_discharge_y[id_cell] += (  edge_bnd_factor_3_F2[edge_idx] * F2_bnd[edge_idx]
                                      + edge_bnd_factor_3_S2[edge_idx] * s2R_bnd[edge_idx]
                                      + edge_bnd_factor_3_F3[edge_idx] * F3_bnd[edge_idx]) * dt;
    }
    //==== Boundary edges
     
    //---- Internal edges
    for (size_t edge_idx = 0; edge_idx < edge_int_nb ; edge_idx++)
    {
        size_t id_celL = edge_int_celL_idx[edge_idx];
        size_t id_celR = edge_int_celR_idx[edge_idx];
        water_depth[id_celL] += edge_int_factor_1_F1_L[edge_idx] * F1_int[edge_idx] * dt;
        water_depth[id_celR] += edge_int_factor_1_F1_R[edge_idx] * F1_int[edge_idx] * dt;
        unit_discharge_x[id_celL] += (  edge_int_factor_2_F2_L[edge_idx] * F2_int[edge_idx]
                                      + edge_int_factor_2_S2_L[edge_idx] * s2L_int[edge_idx]
                                      + edge_int_factor_2_F3_L[edge_idx] * F3_int[edge_idx]) * dt;
        unit_discharge_x[id_celR] += (  edge_int_factor_2_F2_R[edge_idx] * F2_int[edge_idx]
                                      + edge_int_factor_2_S2_R[edge_idx] * s2R_int[edge_idx]
                                      + edge_int_factor_2_F3_R[edge_idx] * F3_int[edge_idx]) * dt;
        unit_discharge_y[id_celL] += (  edge_int_factor_3_F2_L[edge_idx] * F2_int[edge_idx]
                                      + edge_int_factor_3_S2_L[edge_idx] * s2L_int[edge_idx]
                                      + edge_int_factor_3_F3_L[edge_idx] * F3_int[edge_idx]) * dt;
        unit_discharge_y[id_celR] += (  edge_int_factor_3_F2_R[edge_idx] * F2_int[edge_idx]
                                      + edge_int_factor_3_S2_R[edge_idx] * s2R_int[edge_idx]
                                      + edge_int_factor_3_F3_R[edge_idx] * F3_int[edge_idx]) * dt;
    }
    //==== Internal edges
    //== compute balance

    //-- Compute momentum dissipation
    if (model_momentum_dissipation == "active")
    {
        for (size_t idx_cell = 0; idx_cell < cell_nb ; idx_cell++)
        {
            if (water_depth[idx_cell] > water_depth_old[idx_cell] + small_real) 
            //-- Momentum dissipation only if h is increasing
            //-- If mu_1^2 + mu_2^2 = 0 then the flag should be set to 'inactive'
            //-- mu is a dissipation coefficient: mu=0 means that there is no dissipation
            {
                // Get the discharge increment
                double delta_q = unit_discharge_x[idx_cell] - unit_discharge_x_old[idx_cell];
                double delta_r = unit_discharge_y[idx_cell] - unit_discharge_y_old[idx_cell];
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
                unit_discharge_x[idx_cell] -= delta_q_reduction;
                unit_discharge_y[idx_cell] -= delta_r_reduction;
            }
        }
    }
    //== Compute momentum dissipation


    //-- check balance
    //---- Rounding-off error
    for (size_t n = 0; n < cell_nb; ++n)
    {
        if (std::fabs(water_depth[n]) < small_real) water_depth[n] = 0.;
        if (std::fabs(unit_discharge_x[n]) < small_real) unit_discharge_x[n] = 0.;
        if (std::fabs(unit_discharge_y[n]) < small_real) unit_discharge_y[n] = 0.;
    }
    //==== Rounding-off error

    //== check balance
}


void swDip::hypComputeFluxes(void)
/**
 * @brief Computes the fluxes by solving the Riemann problems
 * @return nothing
 * @par Modified attributes:
 *  - all fluxes
 */
{
    //-- Boundary edges
    for (swBoundaryConditionTimeSeries* bc :bc_timeseries_list)
    {
        bc->computeValBC(t);
    }
    for (size_t edge_idx = 0; edge_idx < edge_bnd_nb; edge_idx++)
    {
        double zR = z_riemann_bnd[edge_idx];
        double qR_rotated = q_riemann_bnd[edge_idx];
        double rR_rotated = r_riemann_bnd[edge_idx];

        double flux1 = 0.0;
        double flux2 = 0.0;
        double flux3 = 0.0;
        double source_right = 0.0;

        /** @todo manage rating curves */
        // Get BC type and value
        double bc_val = edge_bnd_bc_pointer[edge_idx]->getBCVal();
        std::string bc_type = edge_bnd_bc_pointer[edge_idx]->getBCType();
        size_t cell_idx = edge_bnd_cell_idx[edge_idx];
        hypBoundaryFlux(boussinesq[cell_idx], phiG_bnd[edge_idx], phiW[cell_idx], cell_z_b[cell_idx],
                        zR, qR_rotated, rR_rotated, 
                        lambda_plus_bnd[edge_idx], bc_type, bc_val,
                        flux1, flux2, flux3, source_right);

        F1_bnd[edge_idx] = flux1;
        F2_bnd[edge_idx] = flux2;
        F3_bnd[edge_idx] = flux3;
        s2R_bnd[edge_idx] = source_right;
    }
    //=== Boundary edges

    //-- Internal edges
    for (size_t edge_idx = 0 ; edge_idx < edge_int_nb ; edge_idx++)
    {
        size_t celL_idx = edge_int_celL_idx[edge_idx];
        size_t celR_idx = edge_int_celR_idx[edge_idx];
        double bL = boussinesq[celL_idx];
        double bR = boussinesq[celR_idx];
        double phiG = phiG_int[edge_idx];
        double phiL = phiW[celL_idx];
        double phiR = phiW[celR_idx];
        double zbL = cell_z_b[edge_int_celL_idx[edge_idx]];
        double zbR = cell_z_b[edge_int_celR_idx[edge_idx]];
        double zL = z_riemann_int[edge_idx][0];
        double zR = z_riemann_int[edge_idx][1];
        double hL = zL - zbL;
        double hR = zR - zbR;
        double qL = q_riemann_int[edge_idx][0];
        double qR = q_riemann_int[edge_idx][1];
        double rL = r_riemann_int[edge_idx][0];
        double rR = r_riemann_int[edge_idx][1];

        // Compute fluxes only if sufficient depth

        if ((hL > h_min) || (hR > h_min))
        {
            double flux1 = 0;
            double flux2 = 0;
            double flux3 = 0;
            double source_left = 0;
            double source_right = 0;

            /** @todo use pointers here */
            std::array<double, 2> lambda_plus_LR = lambda_plus_int[edge_idx];
            std::array<double, 2> lambda_minus_LR = lambda_minus_int[edge_idx];
            if (topo_discretization == "discontinuous")
            {
                hypHllc(bL, bR, phiG, phiL, phiR, zbL, zbR,
                        zL, qL, rL, zR, qR, rR,
                        lambda_minus_LR[0], lambda_minus_LR[1],
                        lambda_plus_LR[0], lambda_plus_LR[1],
                        flux1, flux2, flux3,
                        source_left, source_right);
            }
            else if (topo_discretization == "blending")
            {
                hypHllcBlending(bL, bR, phiG, phiL, phiR, zbL, zbR,
                                zL, qL, rL, zR, qR, rR,
                                lambda_minus_LR[0], lambda_minus_LR[1],
                                lambda_plus_LR[0], lambda_plus_LR[1],
                                flux1, flux2, flux3,
                                source_left, source_right);
            }
            else
            {
                logError() << "'topo_discretization' parameter value not recognized" << topo_discretization;
                logFatal(__FILE__, __func__, __LINE__);
            }

            if (std::isnan(flux1) == true)
            {
                logError() << "DIP - Mass flux is Not a number for internal interface " << edge_idx;
                logInfo() << "DEBUG INTERFACE" << edge_int_idx[edge_idx]
                    << "nx=" << edge_int_norm_x[edge_idx]
                    << "ny=" << edge_int_norm_y[edge_idx];
                logInfo() << "  Left-hand cell" << edge_int_celL_idx[edge_idx]
                    << "zb=" << zbL
                    << "z=" << zL << "q_rotated=" << qL << "r_rotated=" << rL;
                logInfo() << "  Right-hand cell " << edge_int_celR_idx[edge_idx]
                    << "zb=" << zbR
                    << "z=" << zR << "q_rotated=" << qR << "r_rotated=" << rR;
                logInfo() << "  flux1=" << flux1
                    << "flux2=" << flux2
                    << "flux3=" << flux3
                    << "S2L=" << source_left
                    << "S2R=" << source_right;
                logFatal(__FILE__, __func__, __LINE__);
            }

            F1_int[edge_idx] = flux1;
            F2_int[edge_idx] = flux2;
            F3_int[edge_idx] = flux3;
            s2L_int[edge_idx] = source_left;
            s2R_int[edge_idx] = source_right;
        }
    }
    //== Internal edges
}


void swDip::hypComputeLRStates(std::string const& num_scheme)
/**
 * @brief Defines the local Riemann problems for the given numerical scheme
 * @param[in] num_scheme Numerical scheme used (e.g. godunov, muscl-evr)
 * @return nothing
 * @par Modified attributes:
 *  @todo to be checked
 */
{
    if (num_scheme == "godunov")
    {
        hypGodunovComputeLRStates();
    }
    else if (num_scheme == "muscl_evr")
    {
        hypMusclEvrComputeLRStates();
    }
    else
    {
        logError() << "Numerical scheme is not implemented : " << num_scheme;
        logFatal(__FILE__, __func__, __LINE__);
    }
}


double swDip::hypComputeTimeStep(void)
/**
 * @brief Computes the CFL-permissible computational timestep for the DIP model
 * @return dt : computational time step
 */
{
    double dt = dt_max;
    std::vector<double> swept(cell_nb, 0.);

    // --- Boundary edges
    for (size_t edge_idx = 0; edge_idx < edge_bnd_nb; edge_idx++)
    {
        size_t id_cell = edge_bnd_cell_idx[edge_idx];
        double lambda = 0.;
        double h1 = water_depth[id_cell];
        if (h1 > h_min)
        {
            lambda = lambda_max_bnd[edge_idx];
            swept[id_cell] += lambda * edge_bnd_length[edge_idx];
        }
    }
    // === Boundary edges

    // --- Internal edges
    for (size_t edge_idx = 0; edge_idx < edge_int_nb; edge_idx++)
    {
        double lambda = 0.;
        
        // --- Left-hand side
        size_t id_celL = edge_int_celL_idx[edge_idx];
        double h1 = water_depth[id_celL];
        if (h1 > h_min)
        {
            lambda = lambda_max_int[edge_idx];
        }
        swept[id_celL] += lambda * edge_int_length[edge_idx];
        // === Left-hand side

        // --- Right-hand side
        size_t id_celR = edge_int_celR_idx[edge_idx];
        h1 = water_depth[id_celR];
        if (h1 > h_min)
        {
            lambda = lambda_max_int[edge_idx];
        }
        // === Right-hand side
        swept[id_celR] += lambda * edge_int_length[edge_idx];
    }
    // === Internal edges

    // Cells : determine the maximum Dt
    for (size_t cell_idx = 0; cell_idx < cell_nb; cell_idx++)
    {
        double xxx = cfl_max * (cell_area[cell_id[cell_idx]] /
            std::max(swept[cell_id[cell_idx]], small_real));
        if (xxx < dt)
        {
            dt = xxx;
            if (dt < 0.)
            {
                logError() << "Negative Time step in CFL limiting!  deltat :" << dt << "in cell:" << cell_id[cell_idx];
            }
        }
    }

    //-- Map storage-related time step limitation
    if (dt_map > small_real)
    {
        int t_now_over_dt_int = (int)(t / dt_map);
        int t_later_over_dt_int = (int)((t + dt) / dt_map);
        if (t_later_over_dt_int > t_now_over_dt_int) /* Storing is mandatory */
        {
            dt = t_later_over_dt_int * dt_map - t; /* Revise dt so as to reach exactly the nearest storage time */
        }
    }
    if ((t + dt) > t_max)
    {
        dt = t_max - t;
    }
    //== Map storage-related time step limitation

    return dt;
}


void swDip::hypDivCorr(std::string const& div_corr_option, size_t n_it_div_max, double& dt)
/**
 * @brief Corrects fluxes or time step from divergence considerations
 * @param[in] div_corr_option The desired option for divergence correction
 * @param[in,out] dt Computational time step
 * @return nothing
 */
{
    if (div_corr_option == "limit_timestep")
    {
        hypDivCorrTimeStep();
    }
    else if (div_corr_option == "limit_flux")
    {
        hypDivCorrFluxes(n_it_div_max);
    }
    else
    {
        logError() << "Divergence correction option is not handled: " << div_corr_option;
        logFatal(__FILE__, __func__, __LINE__);
    }
}


void swDip::hypDivCorrFluxes(size_t n_it_div_max)
/**
 * @brief compute the divergence correction (with the fluxes limitation) to avoid overdrying
 * @param[in] n_it_div_max : max number of iterations in the process
 * @par Modified attributes :
 *  - all fluxes
 */
{
    bool correct = false;
    size_t nbiter = 0;

    while (correct == false)
    {
        size_t ndry = 0;
        for (size_t cell_idx = 0 ; cell_idx < cell_nb ; cell_idx++)
        {
            //-- For each cell, compute the flux factor that avoids drying and revise fluxes accordingly
            double flux_ratio = 1.0;
            if (water_depth[cell_idx] < 0.0)
            {
                ndry++;
                flux_ratio = water_depth_old[cell_idx] / (water_depth_old[cell_idx] - water_depth[cell_idx]);
                std::vector<size_t> id_edge_bnd = cell_edges_bnd_idx[cell_idx];
                for (size_t edge_idx = 0; edge_idx < id_edge_bnd.size(); edge_idx++)
                {
                    F1_bnd[edge_idx] *= flux_ratio;
                    F2_bnd[edge_idx] *= flux_ratio;
                    F3_bnd[edge_idx] *= flux_ratio;
                    s2R_bnd[edge_idx] *= flux_ratio;
                }
                std::vector<size_t> id_edge_int = cell_edges_int_idx[cell_idx];
                for (size_t edge_idx = 0; edge_idx < id_edge_int.size(); edge_idx++)
                {
                    F1_int[edge_idx] *= flux_ratio;
                    F2_int[edge_idx] *= flux_ratio;
                    F3_int[edge_idx] *= flux_ratio;
                    s2L_int[edge_idx] *= flux_ratio;
                    s2R_int[edge_idx] *= flux_ratio;
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
            water_depth[idx_celL] = water_depth_old[idx_celL] + edge_int_factor_1_F1_L[idx_edge_int] * F1_int[idx_edge_int] * dt;
            water_depth[idx_celR] = water_depth_old[idx_celR] + edge_int_factor_1_F1_R[idx_edge_int] * F1_int[idx_edge_int] * dt;
            unit_discharge_x[idx_celL] = unit_discharge_x_old[idx_celL] + (  edge_int_factor_2_F2_L[idx_edge_int] * F2_int[idx_edge_int]
                                        + edge_int_factor_2_S2_L[idx_edge_int] * s2L_int[idx_edge_int]
                                        + edge_int_factor_2_F3_L[idx_edge_int] * F3_int[idx_edge_int]) * dt;
            unit_discharge_x[idx_celR] = unit_discharge_x_old[idx_celR] + (  edge_int_factor_2_F2_R[idx_edge_int] * F2_int[idx_edge_int]
                                        + edge_int_factor_2_S2_R[idx_edge_int] * s2R_int[idx_edge_int]
                                        + edge_int_factor_2_F3_R[idx_edge_int] * F3_int[idx_edge_int]) * dt;
            unit_discharge_y[idx_celL] = unit_discharge_y_old[idx_celL] + (  edge_int_factor_3_F2_L[idx_edge_int] * F2_int[idx_edge_int]
                                        + edge_int_factor_3_S2_L[idx_edge_int] * s2L_int[idx_edge_int]
                                        + edge_int_factor_3_F3_L[idx_edge_int] * F3_int[idx_edge_int]) * dt;
            unit_discharge_y[idx_celR] = unit_discharge_y_old[idx_celR] +  (  edge_int_factor_3_F2_R[idx_edge_int] * F2_int[idx_edge_int]
                                        + edge_int_factor_3_S2_R[idx_edge_int] * s2R_int[idx_edge_int]
                                        + edge_int_factor_3_F3_R[idx_edge_int] * F3_int[idx_edge_int]) * dt;
        }
        //== Recompute balance
        if ((ndry == 0) || (nbiter > n_it_div_max))
        {
            correct = true;
        }
    }
}


void swDip::hypDivCorrTimeStep(void)
/**
 * @brief Corrects time step to prevent overdrying
 * @par Modified attributes :
 *  - dt
 */
{
    double dt_initial = dt;
    double dt_avoid_drying = dt;
    bool some_cells_will_be_dry = false;
    //-- Check whether some cells will become dry, correct dt as appropriate
    for (size_t cell_idx = 0; cell_idx < cell_nb; cell_idx++)
    {
        size_t id_cell = cell_id[cell_idx];
        if (water_depth[id_cell] <= 0.0)
        {
            double dt_avoid_drying_cell = dt * (water_depth_old[id_cell] - water_depth[id_cell]) / water_depth[id_cell];
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
        for (size_t cell_idx = 0; cell_idx < cell_nb; cell_idx++)
        {
            size_t id_cell = cell_id[cell_idx];
            water_depth[id_cell] = water_depth_old[id_cell] * (1.0 - dt_ratio) + water_depth[id_cell] * dt_ratio;
            unit_discharge_x[id_cell] = unit_discharge_x_old[id_cell] * (1.0 - dt_ratio) + unit_discharge_x[id_cell] * dt_ratio;
            unit_discharge_y[id_cell] = unit_discharge_y_old[id_cell] * (1.0 - dt_ratio) + unit_discharge_y[id_cell] * dt_ratio;
        }
        //== Recompute variables at next time level if needed
    }
}


double swDip::hypHCriticalPhiq(double const& q, double const& phiG, double const& PhiW)
/**
 * @brief swDipHyp::hCritical
 * @param[in] q unit discharge
 * @param[in] phiG interface porosity
 * @param[in] PhiW  cell porosity
 * @return
 */
{
    double h_critical = q * PhiW / std::max(phiG, small_real);
    //-- @todo: remove once proven ok:  h_critical = h_critical * h_critical / m_param->g;
    h_critical = h_critical * h_critical / gravity;
    h_critical = std::exp((1./3.) * std::log(h_critical));
    return h_critical;
}


double swDip::hypHCriticalPhiWq(double const& phiwq, double const& phiG)
/**
 * @brief swDipHyp::hCritical
 * @param[in] phiwq 
 * @param[in] phiG interface porosity
 * @return
 */
{
    double h_critical = phiwq / std::max(phiG, small_real);
    //-- @todo: remove once proven ok:  h_critical = h_critical * h_critical / m_param->g;
    h_critical = h_critical * h_critical / gravity;
    h_critical = std::exp((1./3.) * std::log(h_critical));
    return h_critical;
}


void swDip::hypHllc(double const& bL, double const& bR,
						 double const& phiG, double const& phiL, double const& phiR,
						 double const& zbL, double const& zbR,
                         double const& zL, double const& qL, double const& rL,
                         double const& zR, double const& qR, double const& rR,
                         const double& lambda_minus_L, const double& lambda_minus_R,
                         const double& lambda_plus_L, const double& lambda_plus_R, 
                         double& F1, double& F2, double& F3,
                         double& s2L, double& s2R)
/**
 * @brief computes the fluxes using the classical HLLC formulation in Guinot & Soares-Frazao (2006)
 * @param[in] bL Boussinesq coefficient in the cell on the LHS of the interface
 * @param[in] bR Boussinesq coefficient in the cell on the RHS of the interface
 * @param[in] phiL porosity in the cell on the LHS of the interface
 * @param[in] phiR porosity in the cell on the RHS of the interface
 * @param[in] zbL bottom level in the cell on the LHS of the interface
 * @param[in] zbR bottom level in the cell on the RHS of the interface
 * @param[in] zL  free surface elevation in the cell on the LHS of the interface
 * @param[in] qL  normal unit discharge in the cell on the LHS of the interface
 * @param[in] rL  tangent unit discharge in the cell on the LHS of the interface
 * @param[in] zR  free surface elevation in the cell on the RHS of the interface
 * @param[in] qR  normal unit discharge in the cell on the RHS of the interface
 * @param[in] rR  tangent unit discharge in the cell on the RHS of the interface
 * @param[in] lambda_minus_L lambda- value in the cell on the LHS of the interface
 * @param[in] lambda_minus_R lambda- value in the cell on the RHS of the interface
 * @param[in] lambda_plus_L lambda+ value in the cell on the LHS of the interface
 * @param[in] lambda_plus_R lambda+ value in the cell on the RHS of the interface
 * @param[out] flux1 mass flux in the direction normal to the interface
 * @param[out] flux2 momentum flux in the direction normal to the interface (specific force per unit width)
 * @param[out] flux3 momentum flux in the direction tangent to the interface
 * @param[out] source_leftgeometric source term on the LHS of the interface
 * @param[out] source_rightgeometric source term on the RHS of the interface
 */
{

    double x1L = phiL / std::max(phiG, small_real);
    double x1R = phiR / std::max(phiG, small_real);

    double hL = zL - zbL;
    double hR = zR - zbR;

    //-- Identify flow configuration
    bool flow_in_both_cells = false;
    bool flow_in_lhs_cell_alone = false;
    bool flow_in_rhs_cell_alone = false;
    bool flow_nowhere = false;
    if (((hL <= h_min) && (hR <= h_min)) || ((phiL <= phi_min) && (phiR <= phi_min)))
    {
        // Both cells are currently empty, or cannot store any water
        flow_nowhere = true;
    }
    else
    {
        // Flow assumed in both cells by default, assumption to be checked with the next 2 tests
        flow_in_both_cells = true;
        if ((phiR <= phi_min) || ((phiR > phi_min) && (hR <= h_min) && (zL <= zbR + h_min) ))
        {
            // Water cannot flow into RHS cell
            flow_in_lhs_cell_alone = true;
            flow_in_both_cells = false;
        }
        if ((phiL <= phi_min) || ((phiL > phi_min) && (hL <= h_min) && (zR <= zbL + h_min)))
        {
            // Water cannot flow int LHS cell
            flow_in_rhs_cell_alone = true;
            flow_in_both_cells = false;
        }
    }
    size_t sum_bools = int(flow_in_both_cells) + int(flow_in_lhs_cell_alone) + int(flow_in_rhs_cell_alone) + int(flow_nowhere);
    if (sum_bools != 1)
    {
        logError() << "Impossible flow configuration in mpHypHllcAA ";
        logError() << "                   flow_allowed                                   = " << flow_in_both_cells;
        logError() << "                   flow_in_lhs_cell_alone                         = " << flow_in_lhs_cell_alone;
        logError() << "                   flow_in_rhs_cell_alone                         = " << flow_in_rhs_cell_alone;
        logError() << "                   flow_nowhere                                   = " << flow_nowhere;
        logError() << "                   (hL > small_real) && (hR > small_real)         = " << ((hL > small_real) && (hR > small_real));
        logError() << "                   ((phiL <= small_real) && (phiR <= small_real)) = " << ((phiL <= small_real) && (phiR <= small_real));
        logFatal(__FILE__, __func__, __LINE__);
    }
    //== Identify flow configuration

    if (flow_in_both_cells)
    {
        double uL = qL / std::max(hL, small_real);
        double uR = qR / std::max(hR, small_real);

        double L1 = std::min(lambda_minus_L, std::min(lambda_minus_R, 0.));
        double L2 = std::max(lambda_plus_L, std::max(lambda_plus_R, 0.));

        //Compute F1
        F1 = L2 * phiL * qL - L1 * phiR * qR + L1 * L2 * phiG * (zR - zL);
        F1 = F1 / std::max(L2 - L1, small_real);

        // Compute S2L and S2R
        double deltah = std::max(0., std::min(zL, zbR) - zbL);
        double h1 = std::max(zL - zbR, 0.);
        double s2L_b = - gravity * deltah * 0.5 * (h1 + hL);
        s2L_b = phiG * s2L_b;
        double s2L_phi =  .5 * gravity * hL * hL * ( std::min(0., phiG - phiL));

        deltah = std::max(0., std::min(zR, zbL) - zbR);
        h1 = std::max(zR - zbL, 0.);
        double s2R_b = gravity * deltah * 0.5 * (h1 + hR);
        s2R_b = phiG * s2R_b;
        double s2R_phi = .5 * gravity * hR * hR * ( std::max(0., phiR - phiG));

        double s2 = s2L_b + s2R_b;
        s2L = -L1 * s2 / std::max(L2 - L1, small_real) + s2L_phi;
        s2R = L2 * s2 / std::max(L2 - L1, small_real) + s2R_phi;


        // Compute F2
        double FL = 0.;
        double FR = 0.;
        FL = specificForce(bL, phiG, phiL, hL, uL);
        FR = specificForce(bR, phiG, phiR, hR, uR);
        F2 = L2 * FL - L1 * FR + L1 * L2 * (phiR * qR - phiL * qL);
        F2 = F2 / std::max(L2 - L1, small_real);

        //Compute F3
        double vL = x1L * rL / std::max(hL, small_real);
        double vR = x1R * rR / std::max(hR, small_real);
        if (F1 >= 0)
        {
            F3 = F1 * vL;
        }
        else
        {
            F3 = F1 * vR;
        }
    }
    else if (flow_in_rhs_cell_alone)
    {
        //The cell L is empty and the water level in the cell R
        //is below zbL - filling is impossible
        F1 = 0.;
        F2 = 0.;
        s2L = 0.;
        s2R = phiR * 0.5 * gravity * hR * hR;
        F3 = 0.;    // ! Sandra 24/01/2005
    }
    else if (flow_in_lhs_cell_alone)
    {
        //The cell R is empty and the water level in the cell L
        //is below zbR - filling is impossible
        F1 = 0.;
        F2 = 0.;
        s2L = -phiL * 0.5 * gravity * hL * hL;
        s2R = 0.;
        F3 = 0.;    // ! Sandra 24/01/2005
    }
    else if (flow_nowhere)
    {
        F1 = 0.;
        F2 = 0.;
        s2L = 0.;
        s2R = 0.;
        F3 = 0.;
    }
    else
    {
        logError() << "  HLLC_v2 before Stopping for impossible case ";
        logError() << "            zL zR      " << zL << zR;
        logError() << "            qL qR      " << qL << qR;
        logError() << "            zbL R      " << zbL << zbR;
        logError() << "            hL hR      " << hL << hR;
        logError() << "            hL < small " << (hL < small_real);
        logError() << "            hR < small " << (hR < small_real);
        logError() << "            zL <=zbR+s " << (zL < (zbR + small_real));
        logError() << "            F1 F2      " << F1 << F2;
        logError() << "            sL sR      " << s2L << s2R;
        logError() << "RimSol HLLv2 , Impossible case !";
        logFatal(__FILE__, __func__, __LINE__);
    }
}

void swDip::hypHllcBlending(double const& bL, double const& bR,
						    double const& phiG, double const& phiL, double const&phiR,
						    double const& zbL, double const& zbR,
                            double const& zL, double const& qL, double const& rL,
                            double const& zR, double const& qR, double const& rR,
                            const double& lambda_minus_L, const double& lambda_minus_R,
                            const double& lambda_plus_L, const double& lambda_plus_R, 
                            double& F1, double& F2, double& F3,
                            double& s2L, double& s2R)
/**
 * @brief computes the fluxes using the classical HLLC formulation in Guinot & Soares-Frazao (2006)
 * @param[in] bL : Boussinesq coefficient in the cell on the LHS of the interface
 * @param[in] bR : Boussinesq coefficient in the cell on the RHS of the interface
 * @param[in] phiL : porosity in the cell on the LHS of the interface
 * @param[in] phiR : porosity in the cell on the RHS of the interface
 * @param[in] zbL : bottom level in the cell on the LHS of the interface
 * @param[in] zbR : bottom level in the cell on the RHS of the interface
 * @param[in] zL  : free surface elevation in the cell on the LHS of the interface
 * @param[in] qL  : normal unit discharge in the cell on the LHS of the interface
 * @param[in] rL  : tangent unit discharge in the cell on the LHS of the interface
 * @param[in] zR  : free surface elevation in the cell on the RHS of the interface
 * @param[in] qR  : normal unit discharge in the cell on the RHS of the interface
 * @param[in] rR  : tangent unit discharge in the cell on the RHS of the interface
 * @param[in] lambda_minus_L : lambda- value in the cell on the LHS of the interface
 * @param[in] lambda_minus_R : lambda- value in the cell on the RHS of the interface
 * @param[in] lambda_plus_L : lambda+ value in the cell on the LHS of the interface
 * @param[in] lambda_plus_R : lambda+ value in the cell on the RHS of the interface
 * @param[out] F1 : mass flux in the direction normal to the interface
 * @param[out] F2 : momentum flux in the direction normal to the interface (specific force per unit width)
 * @param[out] F3 : momentum flux in the direction tangent to the interface
 * @param[out] s2L: geometric source term on the LHS of the interface
 * @param[out] s2R: geometric source term on the RHS of the interface
 */
{

    double x1L = phiL / std::max(phiG, small_real);
    double x1R = phiR / std::max(phiG, small_real);

    double hL = zL - zbL;
    double hR = zR - zbR;

    //double zbmax = std::max(zbL, zbR);
    if ((hL > small_real && hR > small_real) ||
        (hL < small_real && zR >= zbL) ||
        (hR < small_real && zL >= zbR))
    {
        //Waves speed
        double uL = qL / std::max(hL, small_real);
        double uR = qR / std::max(hR, small_real);
        double L1 = std::min(lambda_minus_L, std::min(lambda_minus_R, 0.));
        double L2 = std::max(lambda_plus_L, std::max(lambda_plus_R, 0.));

        double alpha_h = 0.0;
		if (std::abs(zbL - zbR) > small_real)
		{
			alpha_h = std::max(0.0 , std::min(1.0 , (zL - zR) / (zbL - zbR) ));
		}

        //Compute F1
        double F1h = L2 * phiL * qL - L1 * phiR * qR + L1 * L2 * phiG * (hR - hL);
        double F1z = L2 * phiL * qL - L1 * phiR * qR + L1 * L2 * phiG * (zR - zL);
        F1 = alpha_h * F1h + (1.0 - alpha_h) * F1z;
        F1 = F1 / std::max(L2 - L1, small_real);

        // Compute S2L and S2R
        double deltah = std::max(0., std::min(zL, zbR) - zbL);
        double h1 = std::max(zL - zbR, 0.);
        double s2L_b_z = - gravity * deltah * 0.5 * (h1 + hL);
        s2L_b_z = phiG * s2L_b_z;
        double s2L_b_h = phiG * hL * std::min(0.0, zbL - zbR);
        double s2L_b = alpha_h * s2L_b_h + (1.0 - alpha_h) * s2L_b_z; 
        double s2L_phi =  .5 * gravity * hL * hL * ( std::min(0., phiG - phiL));

        deltah = std::max(0., std::min(zR, zbL) - zbR);
        h1 = std::max(zR - zbL, 0.);
        double s2R_b_z = gravity * deltah * 0.5 * (h1 + hR);
        s2R_b_z = phiG * s2R_b_z;
        double s2R_b_h = phiG * hR * std::max(0.0, zbL - zbR);
        double s2R_b = alpha_h * s2R_b_h + (1.0 - alpha_h) * s2R_b_z;
        double s2R_phi = .5 * gravity * hR * hR * ( std::max(0., phiR - phiG));

        double s2 = s2L_b + s2R_b;
        s2L = -L1 * s2 / std::max(L2 - L1, small_real) + s2L_phi;
        s2R = L2 * s2 / std::max(L2 - L1, small_real) + s2R_phi;

        // Compute F2
        double FL = 0.;
        double FR = 0.;
        FL = specificForce(bL, phiG, phiL, hL, uL);
        FR = specificForce(bR, phiG, phiR, hR, uR);
        F2 = L2 * FL - L1 * FR + L1 * L2 * (phiR * qR - phiL * qL);
        F2 = F2 / std::max(L2 - L1, small_real);

        //Compute F3
        double vL = x1L * rL / std::max(hL, small_real);
        double vR = x1R * rR / std::max(hR, small_real);
        if (F1 >= 0)
        {
            F3 = F1 * vL;
        }
        else
        {
            F3 = F1 * vR;
        }
        /*
        if (debug_solver == 1) {
            logInfo() << "hllc_v2 : PhiWL R  " << phiL << phiR;
            logInfo() << "          phiG     " << phiG;
            logInfo() << "          hL qL rL " << hL << qL << rL;
            logInfo() << "          hR qR rR " << hR << qR << rR;
            logInfo() << "          Lm Lp    " << L1 << L2;
            logInfo() << "          FL       " << phiL * qL << FL * phiG << F3;
            logInfo() << "          FR       " << phiR * qR << FR * phiG << F3;
            logInfo() << "          F        " << F1 << F2 << F3;
            logInfo() << "          s2L R    " << s2L << s2R;
        }
        */
    }
    else if ((hL <= small_real) && (zR <= zbL + small_real))
    {
        //The cell L is empty and the water level in the cell R
        //is below zbL - filling is impossible
        F1 = 0.;
        F2 = 0.;
        s2L = 0.;
        s2R = phiR * 0.5 * gravity * hR * hR;
        F3 = 0.;    // ! Sandra 24/01/2005
    }
    else if ((hR <= small_real) && (zL <= zbR + small_real))
    {
        //The cell R is empty and the water level in the cell L
        //is below zbR - filling is impossible
        F1 = 0.;
        F2 = 0.;
        s2L = -phiL * 0.5 * gravity * hL * hL;
        s2R = 0.;
        F3 = 0.;    // ! Sandra 24/01/2005
    }
    else
    {
        logError() << "  HLLC_v2 before Stopping for impossible case ";
        logError() << "            zL zR      " << zL << zR;
        logError() << "            qL qR      " << qL << qR;
        logError() << "            zbL R      " << zbL << zbR;
        logError() << "            hL hR      " << hL << hR;
        logError() << "            hL < small " << (hL < small_real);
        logError() << "            hR < small " << (hR < small_real);
        logError() << "            zL <=zbR+s " << (zL < (zbR + small_real));
        logError() << "            F1 F2      " << F1 << F2;
        logError() << "            sL sR      " << s2L << s2R;
        logError() << "RimSol HLLv2 , Impossible case !";
        logFatal(__FILE__, __func__, __LINE__);
    }
}


/** @todo: adapt code */

void swDip::hypHllcImpervious(double const& phiG, double const& phiL, double const&phiR,
						           double const& bL, double const& bR,
						           double const& zbL, double const& zbR,
                                   double const& zL, double const& qL,
                                   double const& zR, double const& qR,
                                   const double& lambda_minus_L, const double& lambda_plus_R, 
                                   double& F1, double& F2, double& F3,
                                   double& s2L, double& s2R)    
/**
 * @brief Compute the fluxes and source terms through an interface when its porosity is null
 * @param[in] celL left-hand cell number
 * @param[in] celR right-hand cell number
 * @param[in] i interface number
 * @param[in] zL surface elevation of the left-hand cell w.r.t. the interface;
 * @param[in] qL x unit discharge of the left-hand cell w.r.t. the interface;
 * @param[in] rL y unit discharge of the left-hand cell w.r.t. the interface;
 * @param[in] zR surface elevation of the right-hand cell w.r.t. the interface;
 * @param[in] qR x unit discharge of the right-hand cell w.r.t. the interface;
 * @param[in] rR y unit discharge of the right-hand cell w.r.t. the interface;
 * @param[in] lambda_minus_L  lambda- value in the cell on the LHS of the interface
 * @param[in] lambda_plus_R  lambda+ value in the cell on the RHS of the interface
 * @param[out] F1 mass flux through the interface (output)
 * @param[out] F2 x momentum through the interface (output)
 * @param[out] F3 y momentum through the interface (output)
 * @param[out] s2L left-distributed source term (bottom and porosity variations) at the interface (output)
 * @param[out] s2R right-distributed source term (bottom and porosity variations) at the interface (output)
 */
{

    double hL = zL - zbL;
    double hR = zR - zbR;

    //Wave speeds
    double uL = qL / std::max(hL, small_real);
    double uR = qR / std::max(hR, small_real);

    double L1 = std::min(0., lambda_minus_L);
    double L2 = std::max(0., lambda_plus_R);

    //Compute F1, F2, F3
    F1 = 0.;
    F2 = 0.;
    F3 = 0.;

    //Compute S2L and S2R
    s2L = (L1 - uL) * qL - 0.5 * gravity * hL * hL;
    s2L = phiL * s2L;
    s2R = (uR - L2) * qR + 0.5 * gravity * hR * hR;
    s2R = phiR * s2R;
}
