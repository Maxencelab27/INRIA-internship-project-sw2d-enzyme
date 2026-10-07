#include <algorithm>
#include <cmath>
#include <vector>

#include "./swSwes.h"
#include "../common/swSource.h"
#include "../common/swBoundaryConditionTimeSeries.h"
#include "../../../common/swBoundaryCondition.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"
#include "../../../common/param.h"


void swSwes::hypBoundaryFlux(double const& zbR, double const& zR, double const& qR, double const& rR,
                             double const& lambda_plus_R, 
                             std::string const& bc_type, double const& bc_val,
                             double& flux1, double& flux2, double& flux3, double& source_right)
/**
 * @brief Computes fluxes at boundary edges for the SWES models
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
    double flux2R = specificForce(uR, hR);
    //== compute characteristic variables

    //-- Routine Core        
    if (bc_type == "fr")
    {
        // Critical section
        if ((bc_val < 0.) | (bc_val > 1.))
        {
            logWarn() << "incorrect value for 'fr'-type boundary condition:" << bc_val;
            logWarn() << "  value has to be between 0 and 1";
            logFatal(__FILE__, __func__, __LINE__);
        }
        double c = (uR - 2. * cR) / (-bc_val - 2.);
        c = std::max(c, 0.);
        double u = -bc_val * c;
        double h = c * c / gravity;
        flux1 = u * h;
        flux2 = flux2R + (flux1 - q1) * L1;
    }
    else if (bc_type == "h")
    {
        // Prescribed h
        double h_bc = std::max(bc_val, 0.);
        if (hR > h_min)
        {
            flux1 = q1 + L1 * (h_bc - hR);
            flux2 = flux2R + (flux1 - q1) * L1;
        }
        else
        {
            double c = std::sqrt(gravity * h_bc);
            flux1 = h_bc * c;
            flux2 = specificForce(c, h_bc);
        }
    }
    else if (bc_type == "z")
    {
        // Prescribed z
        double h_bc = std::max(bc_val - zbR, 0.);
        if (hR > h_min)
        {
            flux1 = q1 + L1 * (h_bc - hR);
            flux2 = flux2R + (flux1 - q1) * L1;
        }
        else
        {
            double c = std::sqrt(gravity * h_bc);
            flux1 = h_bc * c;
            flux2 = specificForce(c, h_bc);
        }
    }
    else if ((bc_type == "q") ||
             (bc_type == "interp") ||
             (bc_type == "powerlaw"))
    {
        // Prescribed q
        double lambda_critic = 2 * std::pow(gravity * std::max(bc_val, 0.), 1. / 3.);
        double lambda = std::max(L1, lambda_critic);

        flux1 = bc_val;
        flux2 = flux2R + (flux1 - q1) * lambda;
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
        logInfo() << "  Stopping...";
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
}


void swSwes::hypCompute(std::string const& div_corr_option, size_t n_it_div_max, std::string const& num_scheme, double& dt)
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

    t_1 = clock();
    if (div_corr_option != "none")
    {
        hypDivCorr(div_corr_option, n_it_div_max, dt);
    }
    total_t_clock_hyperb_divcorr += clock() - t_1;

    t_1 = clock();
    hypComputeBalance(dt);
    total_t_clock_hyperb_balance += clock() - t_1;

    //---- Handle negative water depth
    if (negative_depth_detection == "active")
    {
        // Detect negetive water depths and send error
        for (size_t n = 0; n < cell_nb; ++n)
        {
            if (water_depth[n] < 0)
            {
                logError() << "SWES - Negative water depth" << water_depth[n]
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


void swSwes::hypComputeBalance(double const& dt)
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
    //---- Boundary edges
    for (size_t idx_edge_bnd = 0; idx_edge_bnd < edge_bnd_nb; idx_edge_bnd++)
    {
        size_t idx_cell = edge_bnd_cell_idx[idx_edge_bnd];
        water_depth[idx_cell] += edge_bnd_factor_1_F1[idx_edge_bnd] * F1_bnd[idx_edge_bnd] * dt;
        unit_discharge_x[idx_cell] += (  edge_bnd_factor_2_F2[idx_edge_bnd] * F2_bnd[idx_edge_bnd]
                                      + edge_bnd_factor_2_S2[idx_edge_bnd] * s2R_bnd[idx_edge_bnd]
                                      + edge_bnd_factor_2_F3[idx_edge_bnd] * F3_bnd[idx_edge_bnd]) * dt;
        unit_discharge_y[idx_cell] += (  edge_bnd_factor_3_F2[idx_edge_bnd] * F2_bnd[idx_edge_bnd]
                                      + edge_bnd_factor_3_S2[idx_edge_bnd] * s2R_bnd[idx_edge_bnd]
                                      + edge_bnd_factor_3_F3[idx_edge_bnd] * F3_bnd[idx_edge_bnd]) * dt;
    }
    //==== Boundary edges
     
    //---- Internal edges
    for (size_t idx_edge_int = 0; idx_edge_int < edge_int_nb ; idx_edge_int++)
    {
        size_t idx_celL = edge_int_celL_idx[idx_edge_int];
        size_t idx_celR = edge_int_celR_idx[idx_edge_int];
        water_depth[idx_celL] += edge_int_factor_1_F1_L[idx_edge_int] * F1_int[idx_edge_int] * dt;
        water_depth[idx_celR] += edge_int_factor_1_F1_R[idx_edge_int] * F1_int[idx_edge_int] * dt;
        unit_discharge_x[idx_celL] += (  edge_int_factor_2_F2_L[idx_edge_int] * F2_int[idx_edge_int]
                                      + edge_int_factor_2_S2_L[idx_edge_int] * s2L_int[idx_edge_int]
                                      + edge_int_factor_2_F3_L[idx_edge_int] * F3_int[idx_edge_int]) * dt;
        unit_discharge_x[idx_celR] += (  edge_int_factor_2_F2_R[idx_edge_int] * F2_int[idx_edge_int]
                                      + edge_int_factor_2_S2_R[idx_edge_int] * s2R_int[idx_edge_int]
                                      + edge_int_factor_2_F3_R[idx_edge_int] * F3_int[idx_edge_int]) * dt;
        unit_discharge_y[idx_celL] += (  edge_int_factor_3_F2_L[idx_edge_int] * F2_int[idx_edge_int]
                                      + edge_int_factor_3_S2_L[idx_edge_int] * s2L_int[idx_edge_int]
                                      + edge_int_factor_3_F3_L[idx_edge_int] * F3_int[idx_edge_int]) * dt;
        unit_discharge_y[idx_celR] += (  edge_int_factor_3_F2_R[idx_edge_int] * F2_int[idx_edge_int]
                                      + edge_int_factor_3_S2_R[idx_edge_int] * s2R_int[idx_edge_int]
                                      + edge_int_factor_3_F3_R[idx_edge_int] * F3_int[idx_edge_int]) * dt;
    }
    //==== Internal edges

    //== compute balance

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


void swSwes::hypComputeFluxes()
/**
 * @brief Computes the fluxes by solving the Riemann problems
 * @return nothing
 *  * @par Modified attributes:
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
        double qR_rotated = q_riemann_bnd[idx_edge_bnd];
        double rR_rotated = r_riemann_bnd[idx_edge_bnd];

        double flux1 = 0.0;
        double flux2 = 0.0;
        double flux3 = 0.0;
        double source_right = 0.0;

        // Get BC type and value
        std::string bc_type = edge_bnd_bc_pointer[idx_edge_bnd]->getBCType();
        //-- if rating curve => recompute discharge
        if ((bc_type == "interp") ||
            (bc_type == "powerlaw"))
        {
            edge_bnd_bc_pointer[idx_edge_bnd]->computeValBC(zR);
        }
        //== if rating curve => recompute discharge
        double bc_val = edge_bnd_bc_pointer[idx_edge_bnd]->getBCVal();

        hypBoundaryFlux(cell_z_b[edge_bnd_cell_idx[idx_edge_bnd]], zR, qR_rotated, rR_rotated,
                        lambda_plus_bnd[idx_edge_bnd], bc_type, bc_val,
                        flux1, flux2, flux3, source_right);

        F1_bnd[idx_edge_bnd] = flux1;
        F2_bnd[idx_edge_bnd] = flux2;
        F3_bnd[idx_edge_bnd] = flux3;
        s2R_bnd[idx_edge_bnd] = source_right;
    }
    //=== Boundary edges

    //-- Internal edges
    for (size_t idx_edge_int = 0 ; idx_edge_int < edge_int_nb ; idx_edge_int++)
    {
        double zbL = cell_z_b[edge_int_celL_idx[idx_edge_int]];
        double zbR = cell_z_b[edge_int_celR_idx[idx_edge_int]];
        double zL = z_riemann_int[idx_edge_int][0];
        double zR = z_riemann_int[idx_edge_int][1];
        double hL = zL - zbL;
        double hR = zR - zbR;
        double qL = q_riemann_int[idx_edge_int][0];
        double qR = q_riemann_int[idx_edge_int][1];
        double rL = r_riemann_int[idx_edge_int][0];
        double rR = r_riemann_int[idx_edge_int][1];

        // Compute fluxes only if sufficient depth
        double flux1 = 0;
        double flux2 = 0;
        double flux3 = 0;
        double source_left = 0;
        double source_right = 0;
        if ((hL > h_min) || (hR > h_min))
        {
            std::array<double, 2> lambda_plus_LR = lambda_plus_int[idx_edge_int];
            std::array<double, 2> lambda_minus_LR = lambda_minus_int[idx_edge_int];
            if (topo_discretization == "discontinuous")
            {
                hypHllc(zbL, zbR,
                        zL, qL, rL, zR, qR, rR,
                        lambda_minus_LR[0], lambda_minus_LR[1],
                        lambda_plus_LR[0], lambda_plus_LR[1],
                        flux1, flux2, flux3,
                        source_left, source_right);
            }
            else if (topo_discretization == "blending")
            {
                hypHllcBlending(zbL, zbR,
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
                logError() << "Swes - Mass flux is Not a number for internal interface " << idx_edge_int;
                logInfo() << "DEBUG INTERFACE" << edge_int_idx[idx_edge_int]
                    << "nx=" << edge_int_norm_x[idx_edge_int]
                    << "ny=" << edge_int_norm_y[idx_edge_int];
                logInfo() << "  Left-hand cell" << edge_int_celL_idx[idx_edge_int]
                    << "zb=" << zbL
                    << "z=" << zL << "q_rotated=" << qL << "r_rotated=" << rL;
                logInfo() << "  Right-hand cell " << edge_int_celR_idx[idx_edge_int]
                    << "zb=" << zbR
                    << "z=" << zR << "q_rotated=" << qR << "r_rotated=" << rR;
                logInfo() << "  flux1=" << flux1
                    << "flux2=" << flux2
                    << "flux3=" << flux3
                    << "S2L=" << source_left
                    << "S2R=" << source_right;
                logFatal(__FILE__, __func__, __LINE__);
            }
        }
        F1_int[idx_edge_int] = flux1;
        F2_int[idx_edge_int] = flux2;
        F3_int[idx_edge_int] = flux3;
        s2L_int[idx_edge_int] = source_left;
        s2R_int[idx_edge_int] = source_right;
    }
    //== Internal edges
}


void swSwes::hypComputeLRStates(std::string const& num_scheme)
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


double swSwes::hypComputeTimeStep(void)
/**
 * @brief Computes the CFL-permissible computational timestep for the SWES model
 * @return dt : computational time step
 */
{
    double dt = dt_max;
    std::vector<double> swept(cell_nb, 0.);

    // --- Boundary edges
    for (size_t idx_edge_bnd = 0; idx_edge_bnd < edge_bnd_nb; idx_edge_bnd++)
    {
        size_t idx_cell = edge_bnd_cell_idx[idx_edge_bnd];
        double lambda = 0.;
        double h1 = water_depth[idx_cell];
        if (h1 > h_min)
        {
            lambda = lambda_max_bnd[idx_edge_bnd];
            swept[idx_cell] += lambda * edge_bnd_length[idx_edge_bnd];
        }
    }
    // === Boundary edges

    // --- Internal edges
    for (size_t idx_edge_int = 0; idx_edge_int < edge_int_nb; idx_edge_int++)
    {
        double lambda = 0.;
        
        // --- Left-hand side
        size_t idx_celL = edge_int_celL_idx[idx_edge_int];
        double h1 = water_depth[idx_celL];
        if (h1 > h_min)
        {
            lambda = lambda_max_int[idx_edge_int];
        }
        swept[idx_celL] += lambda * edge_int_length[idx_edge_int];
        // === Left-hand side

        // --- Right-hand side
        size_t idx_celR = edge_int_celR_idx[idx_edge_int];
        h1 = water_depth[idx_celR];
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
                logError() << "Negative Time step in Courant limiting!  deltat :" << dt << "in cell:" << cell_id[idx_cell];
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


void swSwes::hypDivCorr(std::string const& div_corr_option, size_t n_it_div_max, double& dt)
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


void swSwes::hypDivCorrFluxes(size_t n_it_div_max)
/**
 * @brief compute the divergence correction (with the fluxes limitation) to avoid overdrying
 * @param[in] n_it_div_max : max number of iterations in the process
 * @par Modified attributes :
 *  - all fluxes
 */
{
    bool correct = false;
    size_t nbiter = 0;
    double ratio = 1.;

    while (correct == false)
    {
        size_t ndry = 0;
        for (size_t idx_cell = 0 ; idx_cell < cell_nb ; idx_cell++)
        {
            double dh = 0.;
            ///-- Compute the total depth variation over the time step
            for (size_t idx = 0 ; idx < cell_edges_bnd_idx[idx_cell].size() ; idx++)
            {
                size_t idx_edge_bnd = cell_edges_bnd_idx[idx_cell][idx];
                double delta_dh = dt * edge_bnd_factor_1_F1[idx_edge_bnd] * F1_bnd[idx_edge_bnd];
                dh = dh + delta_dh;
                if (std::isnan(dh) == true)
                {
                    logError() << "NAN water depth variation: dh = " << dh;
                    logError() << "  Boundary interface i = " << idx_edge_bnd;
                    logFatal(__FILE__, __func__, __LINE__);
                }
            }
            for (size_t idx = 0 ; idx < cell_edges_int_idx[idx_cell].size() ; idx++)
            {
                size_t idx_edge_int = cell_edges_int_idx[idx_cell][idx];
                double delta_dh = dt * edge_int_factor_1_F1_R[idx_edge_int] * F1_int[idx_edge_int];
                if (idx_cell == edge_int_celL_idx[idx_edge_int])
                {
                    dh = dh - delta_dh;
                }
                else if (idx_cell == edge_int_celR_idx[idx_edge_int])
                {
                    dh = dh + delta_dh;
                }
                else 
                {
                    logError() << "Cell" << idx_cell << "and edge" << idx_edge_int << "do not correspond.";
                    logFatal(__FILE__, __func__, __LINE__);
                }
                if (std::isnan(dh) == true)
                {
                    logError() << "NAN water depth variation: dh = " << dh;
                    logError() << "  Internal interface i = " << idx_edge_int;
                    logFatal(__FILE__, __func__, __LINE__);
                }
            }
            //== Compute the total depth variation over the time step

            ///-- If dh too large, then reduce the fluxes and source terms
            // over all the edges of the current cell
            if ((water_depth[idx_cell] > h_min) &&
                    (dh < - water_depth[idx_cell]))
            {
                ++ndry;
                ratio = - water_depth[idx_cell] / dh;
                for (size_t idx = 0 ; idx < cell_edges_bnd_idx[idx_cell].size() ; idx++)
                {
                    size_t idx_edge_bnd = cell_edges_bnd_idx[idx_cell][idx];
                    /*
                        Correct only boundary edges from which
                        the water is leaving. Inflow boundaries are left unchanged
                    */
                    if (F1_bnd[idx_edge_bnd] < - small_real)
                    {
                        F1_bnd[idx_edge_bnd]  = F1_bnd[idx_edge_bnd] * ratio;
                        F2_bnd[idx_edge_bnd]  = F2_bnd[idx_edge_bnd] * ratio;
                        F3_bnd[idx_edge_bnd]  = F3_bnd[idx_edge_bnd] * ratio;
                        s2R_bnd[idx_edge_bnd] = s2R_bnd[idx_edge_bnd] * ratio;
                    }
                }
                for (size_t idx = 0 ; idx < cell_edges_int_idx[idx_cell].size() ; idx++)
                {
                    size_t idx_edge_int = cell_edges_int_idx[idx_cell][idx];
                    F1_int[idx_edge_int]  = F1_int[idx_edge_int] * ratio;
                    F2_int[idx_edge_int]  = F2_int[idx_edge_int] * ratio;
                    F3_int[idx_edge_int]  = F3_int[idx_edge_int] * ratio;
                    s2L_int[idx_edge_int] = s2L_int[idx_edge_int] * ratio;
                    s2R_int[idx_edge_int] = s2R_int[idx_edge_int] * ratio;
                }
            }
            //== If dh too large, then reduce the fluxes and source terms
        }
        ++nbiter;
        if ((ndry == 0) || (nbiter > n_it_div_max))
        {
            correct = true;
        }
    }
}


void swSwes::hypDivCorrTimeStep(void)
/**
 * @brief Corrects time step to prevent overdrying
 * @par Modified attributes :
 *  - dt
 */
{
    double ratio = 1.;
    for (size_t idx_cell = 0 ; idx_cell < cell_nb ; idx_cell++)
    {
        double dh = 0.;
        // Compute the total depth variation over the time step
        if (water_depth[idx_cell] > h_min)
        {
            for (size_t idx = 0 ; idx < cell_edges_bnd_idx[idx_cell].size() ; idx++)
            {
                size_t idx_edge_bnd = cell_edges_bnd_idx[idx_cell][idx];
                double delta_dh = dt * edge_bnd_factor_1_F1[idx_edge_bnd] * F1_bnd[idx_edge_bnd];
                dh = dh + delta_dh;
                if (std::isnan(dh) == true)
                {
                    logError() << "NAN water depth variation: dh = " << dh;
                    logError() << "  Boundary interface i = " << idx_edge_bnd;
                    logFatal(__FILE__, __func__, __LINE__);
                }
            }
            for (size_t idx = 0 ; idx < cell_edges_int_idx[idx_cell].size() ; idx++)
            {
                size_t idx_edge_int = cell_edges_int_idx[idx_cell][idx];
                double delta_dh = dt * edge_int_factor_1_F1_R[idx_edge_int] * F1_int[idx_edge_int];
                if (idx_cell == edge_int_celL_idx[idx_edge_int])
                {
                    // cell is one left-hand side of an internal edge
                    dh = dh - delta_dh;
                }
                else
                {
                    // cell is one right-hand side of an internal edge
                    dh = dh + delta_dh;
                }
                if (std::isnan(dh) == true)
                {
                    logError() << "NAN water depth variation: dh = " << dh;
                    logError() << "  Internal interface i = " << idx_edge_int;
                    logFatal(__FILE__, __func__, __LINE__);
                }
            }

            //If dh too large, then reduce the timestep
            if (dh < - water_depth[idx_cell])
            {
                /* dh is negative and the leaving volume is bigger than the stock: water_depth[cell->getId()] */
                ratio = std::min(ratio, -water_depth[idx_cell] / dh);
            }
        }
    }

    if (ratio < 0.)
    {
        logError() << "Negative ratio during divergence correction process";
        logError() << "Correction Ratio = " << ratio;
        logFatal(__FILE__, __func__, __LINE__);
    }
    dt = dt * ratio;
}


double swSwes::hypHCritical(const double& q)
/**
 * @brief compute h_critical as un function of discharge
 * @param[in] q: unit discharge
 * @return h_critical
 */
{
    double h_critical = std::pow(q * q / gravity, 1. / 3.);
    return h_critical;
}


void swSwes::hypHllc(double const& zbL, double const& zbR,
    double const& zL, double const& qL, double const& rL,
    double const& zR, double const& qR, double const& rR,
    double const& lambda_minus_L, double const& lambda_minus_R,
    double const& lambda_plus_L, double const& lambda_plus_R, 
    double& flux1, double& flux2, double& flux3,
    double& source_left, double& source_right)
/**
 * @brief computes the fluxes using the classical HLLC formulation in Guinot & Soares-Frazao (2006)
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
 * @param[out] flux1 : mass flux in the direction normal to the interface
 * @param[out] flux2 : momentum flux in the direction normal to the interface (specific force per unit width)
 * @param[out] flux3 : momentum flux in the direction tangent to the interface
 * @param[out] source_left: geometric source term on the LHS of the interface
 * @param[out] source_right: geometric source term on the RHS of the interface
 */
{
    double hL = zL - zbL;
    double hR = zR - zbR;

    if ((hL > small_real && hR > small_real) ||
        (hL < small_real && zR >= zbL) ||
        (hR < small_real && zL >= zbR))
    {
        // Waves speed
        double uL = qL / std::max(hL, small_real);
        double uR = qR / std::max(hR, small_real);
        double L1 = std::min(lambda_minus_L, std::min(lambda_minus_R, 0.));
        double L2 = std::max(lambda_plus_L, std::max(lambda_plus_R, 0.));

        // Compute flux1
        flux1 = (L2 * qL - L1 * qR + L1 * L2 * (zR - zL)) / std::max(L2 - L1, small_real);
        // Compute S2L and S2R
        double fact = 0.5 * (hL + hR);
        double s2 = (hR * hR - hL * hL) * 0.5 + fact * (zL - zR);
        s2 = s2 * gravity / std::max(L2 - L1, small_real);
        source_left = -L1 * s2;
        source_right = L2 * s2;

        // Compute flux2
        double FL = specificForce(uL, hL);
        double FR = specificForce(uR, hR);
        flux2 = (L2 * (FL + source_left) - L1 * (FR - source_right) + L1 * L2 * (qR - qL)) / std::max(L2 - L1, small_real);

        // Compute flux3
        double vL = rL / std::max(hL, small_real);
        double vR = rR / std::max(hR, small_real);
        if (flux1 >= 0)
        {
            flux3 = flux1 * vL;
        }
        else
        {
            flux3 = flux1 * vR;
        }
    }
    else if ((hL <= small_real) && (zR <= zbL + small_real))
    {
        // The cell L is empty and the water level in the cell R
        //  is below zbL - filling is impossible
        double uR = qR / std::max(hR, small_real);
        double cR = std::sqrt(gravity * hR);
        double c = cR - 0.5 * uR;
        flux1 = 0.;
        flux2 = 0.;
        source_left = 0.;
        source_right = 0.5 * std::pow(c, 4) / std::max(gravity, small_real);
        flux3 = 0.;    // ! Sandra 24/01/2005
    }
    else if ((hR <= small_real) && (zL <= zbR + small_real))
    {
        // The cell R is empty and the water level in the cell L
        //  is below zbR - filling is impossible
        double uL = qL / std::max(hL, small_real);
        double cL = std::sqrt(gravity * hL);
        double c = 0.5 * uL + cL;
        flux1 = 0.;
        flux2 = 0.;
        source_left = -0.5 * std::pow(c, 4) / std::max(gravity, small_real);
        source_right = 0.;
        flux3 = 0.;    // ! Sandra 24/01/2005
    }
    else
    {
        logError() << "SWES model HLLC solver before Stopping for impossible case ";
        logError() << "            zL zR      " << zL << zR;
        logError() << "            qL qR      " << qL << qR;
        logError() << "            zbL R      " << zbL << zbR;
        logError() << "            hL hR      " << hL << hR;
        logError() << "            hL < small " << (hL < small_real);
        logError() << "            hR < small " << (hR < small_real);
        logError() << "            zL <=zbR+s " << (zL < (zbR + small_real));
        logError() << "            flux1 flux2      " << flux1 << flux2;
        logError() << "            sL sR      " << source_left << source_right;
        logError() << "SWES model HLLC solver, Impossible case !";
        logFatal(__FILE__, __func__, __LINE__);
    }
}


void swSwes::hypHllcBlending(double const& zbL, double const& zbR,
    double const& zL, double const& qL, double const& rL,
    double const& zR, double const& qR, double const& rR,
    double const& lambda_minus_L, double const& lambda_minus_R,
    double const& lambda_plus_L, double const& lambda_plus_R, 
    double& flux1, double& flux2, double& flux3,
    double& source_left, double& source_right)
/**
 * @brief computes the fluxes using the classical HLLC formulation in Guinot & Soares-Frazao (2006)
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
 * @param[out] flux1 : mass flux in the direction normal to the interface
 * @param[out] flux2 : momentum flux in the direction normal to the interface (specific force per unit width)
 * @param[out] flux3 : momentum flux in the direction tangent to the interface
 * @param[out] source_left: geometric source term on the LHS of the interface
 * @param[out] source_right: geometric source term on the RHS of the interface
 */
{
    double hL = zL - zbL;
    double hR = zR - zbR;

    if ((hL > small_real && hR > small_real) ||
        (hL < small_real && zR >= zbL) ||
        (hR < small_real && zL >= zbR))
    {
        // Waves speed
        double uL = qL / std::max(hL, small_real);
        double uR = qR / std::max(hR, small_real);
        double L1 = std::min(lambda_minus_L, std::min(lambda_minus_R, 0.));
        double L2 = std::max(lambda_plus_L, std::max(lambda_plus_R, 0.));
        double alpha_h = 0.0;
        if (std::abs(zbL - zbR) > small_real)
        {
            alpha_h = std::max(0.0, std::min(1.0, (zL - zR) / (zbL - zbR)));
        }

        // Compute flux1
        double flux1h = (L2 * qL - L1 * qR + L1 * L2 * (hR - hL)) / std::max(L2 - L1, small_real);
        double flux1z = (L2 * qL - L1 * qR + L1 * L2 * (zR - zL)) / std::max(L2 - L1, small_real);
        flux1 = alpha_h * flux1h + (1.0 - alpha_h) * flux1z;

        // Compute S2L and S2R
        double fact = 0.5 * (hL + hR);
        double s2z = (hR * hR - hL * hL) * 0.5 + fact * (zL - zR);
        double s2h = 0.0;
        if (zbL > zbR)
        {
            s2h = hR * (zbL - zbR);
        }
        else
        {
            s2h = hL * (zbL - zbR);
        }
        double s2 = alpha_h * s2h + (1.0 - alpha_h) * s2z;
        s2 = s2 * gravity / std::max(L2 - L1, small_real);
        source_left = -L1 * s2;
        source_right = L2 * s2;

        // Compute flux2
        double FL = specificForce(uL, hL);
        double FR = specificForce(uR, hR);
        flux2 = (L2 * (FL + source_left) - L1 * (FR - source_right) + L1 * L2 * (qR - qL)) / std::max(L2 - L1, small_real);

        // Compute flux3
        double vL = rL / std::max(hL, small_real);
        double vR = rR / std::max(hR, small_real);
        if (flux1 >= 0)
        {
            flux3 = flux1 * vL;
        }
        else
        {
            flux3 = flux1 * vR;
        }
    }
    else if ((hL <= small_real) && (zR <= zbL + small_real))
    {
        // The cell L is empty and the water level in the cell R
        //  is below zbL - filling is impossible
        double uR = qR / std::max(hR, small_real);
        double cR = std::sqrt(gravity * hR);
        double c = cR - 0.5 * uR;
        flux1 = 0.;
        flux2 = 0.;
        source_left = 0.;
        source_right = 0.5 * std::pow(c, 4) / std::max(gravity, small_real);
        flux3 = 0.;    // ! Sandra 24/01/2005
    }
    else if ((hR <= small_real) && (zL <= zbR + small_real))
    {
        // The cell R is empty and the water level in the cell L
        //  is below zbR - filling is impossible
        double uL = qL / std::max(hL, small_real);
        double cL = std::sqrt(gravity * hL);
        double c = 0.5 * uL + cL;
        flux1 = 0.;
        flux2 = 0.;
        source_left = -0.5 * std::pow(c, 4) / std::max(gravity, small_real);
        source_right = 0.;
        flux3 = 0.;    // ! Sandra 24/01/2005
    }
    else
    {
        logError() << "SWES model HLLC solver before Stopping for impossible case ";
        logError() << "            zL zR      " << zL << zR;
        logError() << "            qL qR      " << qL << qR;
        logError() << "            zbL R      " << zbL << zbR;
        logError() << "            hL hR      " << hL << hR;
        logError() << "            hL < small " << (hL < small_real);
        logError() << "            hR < small " << (hR < small_real);
        logError() << "            zL <=zbR+s " << (zL < (zbR + small_real));
        logError() << "            flux1 flux2      " << flux1 << flux2;
        logError() << "            sL sR      " << source_left << source_right;
        logError() << "SWES model HLLC solver, Impossible case !";
        logFatal(__FILE__, __func__, __LINE__);
    }
}

void swSwes::hypHllcImpervious(double const& zbL, double const& zbR,
    double const& zL, double const& qL,
    double const& zR, double const& qR,
    double const& lambda_minus_L, double const& lambda_plus_R, 
    double& flux1, double& flux2, double& flux3,
    double& source_left, double& source_right)
/**
 * @brief Computes the fluxes and source terms through an interface when its porosity is null
 * @param[in] zbL : bottom level in the cell on the LHS of the interface
 * @param[in] zbR : bottom level in the cell on the RHS of the interface
 * @param[in] zL  : free surface elevation in the cell on the LHS of the interface
 * @param[in] qL  : normal unit discharge in the cell on the LHS of the interface
 * @param[in] zR  : free surface elevation in the cell on the RHS of the interface
 * @param[in] qR  : normal unit discharge in the cell on the RHS of the interface
 * @param[in] lambda_minus_L : lambda- value in the cell on the LHS of the interface
 * @param[in] lambda_plus_R : lambda+ value in the cell on the RHS of the interface
 * @param[out] flux1 : mass flux in the direction normal to the interface
 * @param[out] flux2 : momentum flux in the direction normal to the interface (specific force per unit width)
 * @param[out] flux3 : momentum flux in the direction tangent to the interface
 * @param[out] source_left: geometric source term on the LHS of the interface
 * @param[out] source_right: geometric source term on the RHS of the interface
 */
{
    double hL = zL - zbL;
    double hR = zR - zbR;

    // Wave speeds
    double uL = qL / std::max(hL, small_real);
    double uR = qR / std::max(hR, small_real);

    // Compute flux1, flux2, flux3
    flux1 = 0.;
    flux2 = 0.;
    flux3 = 0.;

    // Compute S2L and S2R
    source_left = (lambda_minus_L - uL) * qL - 0.5 * gravity * hL * hL;
    source_right = (uR - lambda_plus_R) * qR + 0.5 * gravity * hR * hR;
}
