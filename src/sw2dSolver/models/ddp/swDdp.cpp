#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstring>
#include <string>
#include <vector>

#include "swDdp.h"
#include "../../../common/param.h"
#include "../../../common/swFile.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"

swDdp::swDdp(void) 
/**
 * @brief Constructor for the swSp object
 */
{
    // Empty on purpose
}


swDdp::~swDdp(void)
/**
 * @brief Destructor for the swSp object
 */
{
    // Empty on purpose
}


void swDdp::checkParameterConsistency(void)
/**
 * @brief checks consistency for all solver parameters
 * @return nothing, but fails if any problem is found
 */
{
    //-- Boussinesq
    if (boussinesq.size() != cell_nb)
    {
        logError() << "Boussinesq parameter contains " << boussinesq.size() << " records, while cell number is " << cell_nb;
        logFatal(__FILE__, __func__, __LINE__);
    }
    else
    {
        for (size_t idx_cell = 0 ; idx_cell < cell_nb ; idx_cell++)
        {
            if (boussinesq[idx_cell] < 1.0)
            {
                logError() << "Cell =" << idx_cell << " - boussinesq coeff. should be >= 1 but is equal to " << boussinesq[idx_cell];
                logFatal(__FILE__, __func__, __LINE__);
            }
        }
    }
    //== Boussinesq

    //-- Friction
    if (model_friction != "none")
    {
        if (friction_coeff_direction_1.size() != cell_nb)
        {
            logError() << "Vector size error: friction_coeff_1.size(), cell_nb = " << friction_coeff_direction_1.size() << cell_nb;
            logFatal(__FILE__, __func__, __LINE__);
        }
        else
        {
            for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
            {
                if (friction_coeff_direction_1[idx_cell] < 0.0)
                {
                    logError() << "Friction_coefficient_1 is negative - cell, value = " << idx_cell << friction_coeff_direction_1[idx_cell];
                }
            }
        }
        if (friction_coeff_direction_2.size() != cell_nb)
        {
            logError() << "Vector size error: friction_coeff_2.size(), cell_nb = " << friction_coeff_direction_2.size() << cell_nb;
            logFatal(__FILE__, __func__, __LINE__);
        }
        else
        {
            for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
            {
                if (friction_coeff_direction_2[idx_cell] < 0.0)
                {
                    logError() << "Friction_coefficient_2 is negative - cell, value = " << idx_cell << friction_coeff_direction_2[idx_cell];
                }
            }
        }
        if (friction_alpha.size() != cell_nb)
        {
            logError() << "Vector size error: friction_alpha.size(), cell_nb = " << friction_alpha.size() << cell_nb;
            logFatal(__FILE__, __func__, __LINE__);
        }
    }
    //== Friction

    //-- Infiltration
    if (model_infiltration != "none")
    {
        if ((model_infiltration == "h_proportional") || (model_infiltration == "horton"))
        {
            if (infiltration_param.size() != cell_nb)
            {
                logError() << "Infiltration model " << model_infiltration << "Vector size error : infiltratation _param.size(), cell_nb = " << infiltration_param.size() << cell_nb;
                logFatal(__FILE__, __func__, __LINE__);
            }

            size_t error_number = 0;
            std::vector<size_t> cell_error_idx = {};
            for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
            {
                if (infiltration_param[idx_cell] < 0.)
                {
                    error_number++;
                    cell_error_idx.push_back(idx_cell);
                }
            }
            if (error_number != 0)
            {
                logError() << "Infiltration model = " << model_infiltration << " parameter is negative in cells = " << cell_error_idx;
                logFatal(__FILE__, __func__, __LINE__);
            }
        }

    }
    //== Infiltration

    //-- Porosity
    //---- 1/3: Cell porosity
    if (phiW_tab.size() != cell_nb)
    {
        logError() << "Cell porosity parameter contains " << phiW_tab.size() << " records, while cell number is " << cell_nb;
        logFatal(__FILE__, __func__, __LINE__);
    }
    else
    {
        for (size_t idx_cell = 0 ; idx_cell < cell_nb ; idx_cell++)
        {
            for (size_t tab_level = 0 ; tab_level < phiW_tab[idx_cell].size() ; tab_level++)
            {
                if ((phiW_tab[idx_cell][tab_level] < 0.0) || (phiW_tab[idx_cell][tab_level] > 1.0))
                {
                    logError() << "Cell =" << idx_cell << " - porosity should belong to [0,1] but is equal to " << phiW_tab[idx_cell][tab_level];
                    logFatal(__FILE__, __func__, __LINE__);
                }
            }
        }
    }
    //==== 1/3: Cell porosity
    //---- 2/3: Edge_bnd porosity
    if (phiG_bnd_tab.size() != edge_bnd_nb)
    {
        logError() << "Edge_bnd porosity parameter contains " << phiG_bnd_tab.size() << " records, while boundary edge number is " << edge_bnd_nb;
        logFatal(__FILE__, __func__, __LINE__);
    }
    else
    {
        for (size_t idx_edge_bnd = 0 ; idx_edge_bnd < edge_bnd_nb ; idx_edge_bnd++)
        {
            for (size_t tab_level = 0 ; tab_level < phiG_bnd_tab[idx_edge_bnd].size() ; tab_level++)
            {
                if ((phiG_bnd_tab[idx_edge_bnd][tab_level] < 0.0) || (phiG_bnd_tab[idx_edge_bnd][tab_level] > 1.0))
                {
                    logError() << "Boundary edge =" << idx_edge_bnd << " - porosity should belong to [0,1] but is equal to " << phiG_bnd_tab[idx_edge_bnd][tab_level];
                    logFatal(__FILE__, __func__, __LINE__);
                }
                if (phiG_bnd_tab[idx_edge_bnd][tab_level] > phiW_tab[edge_bnd_cell_idx[idx_edge_bnd]][tab_level]) 
                {
                    logError() << "Boundary edge =" << idx_edge_bnd << " - level = " << tab_level << " - porosity = " << phiG_bnd_tab[idx_edge_bnd][tab_level] << "larger than cell porosity " <<  phiW_tab[edge_bnd_cell_idx[idx_edge_bnd]][tab_level];
                    logFatal(__FILE__, __func__, __LINE__);
                }
            }
        }
    }
    //==== 2/3: Edge_bnd porosity
    //---- 3/3: Edge_int porosity
    if (phiG_int_tab.size() != edge_int_nb)
    {
        logError() << "Edge_bnd porosity parameter contains " << phiG_int_tab.size() << " records, while boundary edge number is " << edge_int_nb;
        logFatal(__FILE__, __func__, __LINE__);
    }
    else
    {
        for (size_t idx_edge_int = 0 ; idx_edge_int < edge_int_nb ; idx_edge_int++)
        {
            for (size_t tab_level = 0 ; tab_level < phiG_int_tab[idx_edge_int].size() ; tab_level++)
            {
                if ((phiG_int_tab[idx_edge_int][tab_level] < 0.0) || (phiG_int_tab[idx_edge_int][tab_level] > 1.0))
                {
                    logError() << "Internal edge =" << idx_edge_int << " - porosity should belong to [0,1] but is equal to " << phiG_int_tab[idx_edge_int][tab_level];
                    logFatal(__FILE__, __func__, __LINE__);
                }
                if (phiG_int_tab[idx_edge_int][tab_level] > phiW_tab[edge_int_celL_idx[idx_edge_int]][tab_level]) 
                {
                    logError() << "Internal edge =" << idx_edge_int << " - porosity = " << phiG_int_tab[idx_edge_int][tab_level] << "larger than LHS cell porosity " <<  phiW_tab[edge_int_celL_idx[idx_edge_int]][tab_level];
                    logFatal(__FILE__, __func__, __LINE__);
                }
                if (phiG_int_tab[idx_edge_int][tab_level] > phiW_tab[edge_int_celR_idx[idx_edge_int]][tab_level]) 
                {
                    logError() << "Internal edge =" << idx_edge_int << " - porosity = " << phiG_int_tab[idx_edge_int][tab_level] << "larger than RHS cell porosity " <<  phiW_tab[edge_int_celR_idx[idx_edge_int]][tab_level];
                    logFatal(__FILE__, __func__, __LINE__);
                }
            }
        }
    }
    //==== 3/3: Edge_int porosity
    //== Porosity


    //-- Runoff
    if (model_runoff != "none")
    {
        //-- Check precipitation station codes + time series
        //---- Check station codes
        if (precipitation_station_code.size() != cell_nb)
        {
            logError() << "Vector size error: precipitation_station_code.size(), cell_nb = " << precipitation_station_code.size() << cell_nb;
            logFatal(__FILE__, __func__, __LINE__);
        }
        else
        {
            size_t error_number = 0;
            std::vector<size_t> cell_error_idx = {};
            for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
            {
                if (precipitation_station_code[idx_cell] < 1)
                {
                    error_number++;
                    cell_error_idx.push_back(idx_cell);
                }
            }
            if (error_number != 0)
            {
                logError() << "Precipitation station code is zero or smaller in cells = " << cell_error_idx;
                logFatal(__FILE__, __func__, __LINE__);
            }
        }
        //==== Check station codes

        //---- Check time series
        for (size_t idx_time = 1; idx_time < precipitation_times_vector.size(); idx_time++)
        {
            if (precipitation_times_vector[idx_time - 1] > precipitation_times_vector[idx_time])
            {
                logError() << "Precipitation time series - non-increasing times for records = " << idx_time - 1 << idx_time;
                logFatal(__FILE__, __func__, __LINE__);
            }
        }
        //==== Check time series
        //== Check precipitation station codes + time series

        //-- Check parameter values
        if (model_runoff == "infiltration_coefficient" || model_runoff == "runoff_coefficient")
        {
            size_t error_number = 0;
            std::vector<size_t> cell_error_idx = {};
            for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
            {
                if (runoff_param[idx_cell] < 0.0 || runoff_param[idx_cell] > 1.0)
                {
                    error_number++;
                    cell_error_idx.push_back(idx_cell);
                }
            }
            if (error_number != 0)
            {
                logError() << "Infiltrataion / runoff coefficient < 0 or > 1 in cells = " << cell_error_idx;
                logFatal(__FILE__, __func__, __LINE__);
            }

        }
        else if (model_runoff == "scs")
        {
            size_t error_number_ia = 0;
            size_t error_number_s = 0;
            std::vector<size_t> cell_error_ia_idx = {};
            std::vector<size_t> cell_error_s_idx = {};
            for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
            {
                if (scs_i_a[idx_cell] < 0.0)
                {
                    error_number_ia++;
                    cell_error_ia_idx.push_back(idx_cell);
                }
                if (scs_s[idx_cell] < 0.0)
                {
                    error_number_s++;
                    cell_error_s_idx.push_back(idx_cell);
                }
            }
            if (error_number_ia != 0)
            {
                logError() << "SCS - Ia negative in cells = " << cell_error_ia_idx;
            }
            if (error_number_s != 0)
            {
                logError() << "SCS - S negative in cells = " << cell_error_s_idx;
            }
            if ((error_number_ia) != 0 || (error_number_s != 0))
            {
                logFatal(__FILE__, __func__, __LINE__);
            }

        }
        //== Check parameter values


    }
    //== Runoff

    //-- Singular head loss
    if (model_singular_head_loss != "none")
    {
        if (singular_head_loss_coeff_direction_1.size() != cell_nb)
        {
            logError() << "Vector size error: singular_head_loss_coeff_1.size(), cell_nb = " << singular_head_loss_coeff_direction_1.size() << cell_nb;
            logFatal(__FILE__, __func__, __LINE__);
        }
        else
        {
            for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
            {
                if (singular_head_loss_coeff_direction_1[idx_cell] < 0.0)
                {
                    logError() << "singular_head_loss_coefficient_1 is negative - cell, value = " << idx_cell << singular_head_loss_coeff_direction_1[idx_cell];
                }
            }
        }
        if (singular_head_loss_coeff_direction_2.size() != cell_nb)
        {
            logError() << "Vector size error: singular_head_loss_coeff_2.size(), cell_nb = " << singular_head_loss_coeff_direction_2.size() << cell_nb;
            logFatal(__FILE__, __func__, __LINE__);
        }
        else
        {
            for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
            {
                if (singular_head_loss_coeff_direction_2[idx_cell] < 0.0)
                {
                    logError() << "singular_head_loss_coefficient_2 is negative - cell, value = " << idx_cell << singular_head_loss_coeff_direction_2[idx_cell];
                }
            }
        }
        if (singular_head_loss_alpha.size() != cell_nb)
        {
            logError() << "Vector size error: singular_head_loss_alpha.size(), cell_nb = " << singular_head_loss_alpha.size() << cell_nb;
            logFatal(__FILE__, __func__, __LINE__);
        }
    }
    //== Singular head loss
    
    //-- Wind
    if (model_wind != "none")
    {
        //-- Check station codes
        if (wind_station_code.size() != cell_nb)
        {
            logError() << "Vector size error: wind_station_code.size(), cell_nb = " << wind_station_code.size() << cell_nb;
            logFatal(__FILE__, __func__, __LINE__);
        }
        else
        {
            size_t error_number = 0;
            std::vector<size_t> cell_error_idx = {};
            for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
            {
                if (wind_station_code[idx_cell] < 1)
                {
                    error_number++ ;
                    cell_error_idx.push_back(idx_cell);
                }
            }
            if (error_number != 0)
            {
                logError() << "Wind station code is zero or smaller in cells = " << cell_error_idx;
                logFatal(__FILE__, __func__, __LINE__);
            }
        }
        //== Check station codes

        //-- Check time series
        for (size_t idx_time = 1; idx_time < wind_times_vector.size(); idx_time++)
        {
            if (wind_times_vector[idx_time - 1] > wind_times_vector[idx_time])
            {
                logError() << "Wind time series - non-increasing times for records = " << idx_time - 1 << idx_time;
                logFatal(__FILE__, __func__, __LINE__);
            }
        }
        //== Check time series

    }
    //== Check Wind
    //== Write Model logs
}


void swDdp::computeEdgeFactors(void)
/**
 * @brief Computes edge factors for each of the 3 state equations
 * @par Modified attributes:
 * - edge_bnd_factor_1_F1
 * - edge_bnd_factor_2_F2
 * - edge_bnd_factor_2_F3
 * - edge_bnd_factor_2_S2
 * - edge_bnd_factor_3_F2
 * - edge_bnd_factor_3_F3
 * - edge_bnd_factor_3_S2
 * - edge_int_factor_1_F1_L/R
 * - edge_int_factor_2_F2_L/R
 * - edge_int_factor_2_F3_L/R
 * - edge_int_factor_2_S2_L/R
 * - edge_int_factor_3_F2_L/R
 * - edge_int_factor_3_F3_L/R
 * - edge_int_factor_3_S2_L/R
 */
{
    ///--  Boundary edge factors
    edge_bnd_factor_1_F1.resize(edge_bnd_nb);
    edge_bnd_factor_2_F2.resize(edge_bnd_nb);
    edge_bnd_factor_2_F3.resize(edge_bnd_nb);
    edge_bnd_factor_2_S2.resize(edge_bnd_nb);
    edge_bnd_factor_3_F2.resize(edge_bnd_nb);
    edge_bnd_factor_3_F3.resize(edge_bnd_nb);
    edge_bnd_factor_3_S2.resize(edge_bnd_nb);
    for (size_t idx_edge_bnd = 0; idx_edge_bnd < edge_bnd_nb; idx_edge_bnd++)
    {
        double w = edge_bnd_length[idx_edge_bnd];
        size_t cell_id = edge_bnd_cell_idx[idx_edge_bnd];
        double area = cell_area[cell_id];
        double nx = edge_bnd_norm_x[idx_edge_bnd];
        double ny = edge_bnd_norm_y[idx_edge_bnd];
        double coeff = w / area;

        edge_bnd_factor_1_F1[idx_edge_bnd] = + coeff;
        edge_bnd_factor_2_F2[idx_edge_bnd] = + coeff * nx;
        edge_bnd_factor_2_F3[idx_edge_bnd] = - coeff * ny;
        edge_bnd_factor_2_S2[idx_edge_bnd] = + coeff * nx;
        edge_bnd_factor_3_F2[idx_edge_bnd] = + coeff * ny;
        edge_bnd_factor_3_F3[idx_edge_bnd] = + coeff * nx;
        edge_bnd_factor_3_S2[idx_edge_bnd] = + coeff * ny;
    }
    //--   Boundary edge factors

    ///--  Internal edge factors
    edge_int_factor_1_F1_L.resize(edge_int_nb);
    edge_int_factor_1_F1_R.resize(edge_int_nb);
    edge_int_factor_2_F2_L.resize(edge_int_nb);
    edge_int_factor_2_F2_R.resize(edge_int_nb);
    edge_int_factor_2_F3_L.resize(edge_int_nb);
    edge_int_factor_2_F3_R.resize(edge_int_nb);
    edge_int_factor_2_S2_L.resize(edge_int_nb);
    edge_int_factor_2_S2_R.resize(edge_int_nb);
    edge_int_factor_3_F2_L.resize(edge_int_nb);
    edge_int_factor_3_F2_R.resize(edge_int_nb);
    edge_int_factor_3_F3_L.resize(edge_int_nb);
    edge_int_factor_3_F3_R.resize(edge_int_nb);
    edge_int_factor_3_S2_L.resize(edge_int_nb);
    edge_int_factor_3_S2_R.resize(edge_int_nb);
    for (size_t idx_edge_int = 0; idx_edge_int < edge_int_nb; idx_edge_int++)
    {
        double w = edge_int_length[idx_edge_int];
        size_t celL_id = edge_int_celL_idx[idx_edge_int];
        size_t celR_id = edge_int_celR_idx[idx_edge_int];
        double area_L = cell_area[celL_id];
        double area_R = cell_area[celR_id];
        double nx = edge_int_norm_x[idx_edge_int];
        double ny = edge_int_norm_y[idx_edge_int];
        double coeff_L = w / area_L;
        double coeff_R = w / area_R;

        edge_int_factor_1_F1_L[idx_edge_int] = - coeff_L;
        edge_int_factor_1_F1_R[idx_edge_int] = + coeff_R;
        edge_int_factor_2_F2_L[idx_edge_int] = - coeff_L * nx;
        edge_int_factor_2_F2_R[idx_edge_int] = + coeff_R * nx;
        edge_int_factor_2_F3_L[idx_edge_int] = + coeff_L * ny;
        edge_int_factor_2_F3_R[idx_edge_int] = - coeff_R * ny;
        edge_int_factor_2_S2_L[idx_edge_int] = + coeff_L * nx;
        edge_int_factor_2_S2_R[idx_edge_int] = + coeff_R * nx;
        edge_int_factor_3_F2_L[idx_edge_int] = - coeff_L * ny;
        edge_int_factor_3_F2_R[idx_edge_int] = + coeff_R * ny;
        edge_int_factor_3_F3_L[idx_edge_int] = - coeff_L * nx;
        edge_int_factor_3_F3_R[idx_edge_int] = + coeff_R * nx;
        edge_int_factor_3_S2_L[idx_edge_int] = + coeff_L * ny;
        edge_int_factor_3_S2_R[idx_edge_int] = + coeff_R * ny;
    }
    //==   Internal edge factors
}


bool swDdp::computePressureForceOmegaOntoGammaOfZ(const size_t& idx_edge_int, const double& zL, const double& zR, double& force_gamma_L, double& force_gamma_R, double& force_bottom_L, double& force_bottom_R)
/**
 * @brief Compute the pressure force against PhiW from z
 * @param idx_edge_int: inernal edge index
 * @param zL: free surface elevation
 * @param[out] force_gamma_L: left pressure force (water-water)
 * @param[out] force_gamma_R: right pressure force (water-water)
 * @param[out] force_bottom_L: wall reaction on the left-hand side
 * @param[out] force_bottom_R: wall reaction on the right-hand side
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 */
{
    size_t idx_celL = edge_int_celL_idx[idx_edge_int];
    //-- Pressure forces on the water-water boundary
    double force_omega_L = pressureForceOmegaOfZ(idx_celL, zL);
    //== Pressure forces on the water-water boundary
    //-- Wall reaction on the left-hand side
    force_gamma_L = pressureForceGammaOfZ(idx_edge_int, zL);    
    force_bottom_L = std::max(0., force_omega_L - force_gamma_L);
    //== Wall reaction on the left-hand side

    size_t idx_celR = edge_int_celR_idx[idx_edge_int];
    //-- Pressure forces on the water-water boundary
    double force_omega_R = pressureForceOmegaOfZ(idx_celR, zR);
    //== Pressure forces on the water-water boundary
    //-- Wall reaction on the right-hand side
    force_gamma_R = pressureForceGammaOfZ(idx_edge_int, zR);    
    force_bottom_R = std::max(0., force_omega_R - force_gamma_R);
    //== Wall reaction on the right-hand side

    return EXIT_SUCCESS;
}


/** @todo: adapt code */
bool swDdp::correctVelocity(void)
/**
 * @brief Eliminates abnormally large velocities
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 * - unit_discharge_x
 * - unit_discharge_y
 */
{
    //---- Eliminate abnormally large velocities

    /// @todo Implement in a proper way once all parameters are known
    std::vector<size_t> limited_cell_idx;
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        double theta = thetaW[idx_cell];
        if (theta > h_min)
        {
            double q = thetaW_velocity_x[idx_cell];
            double r = thetaW_velocity_y[idx_cell];
            double q_norm = std::sqrt(q * q + r * r);
            double c = std::sqrt(gravity * theta);
            double Fr = q_norm / c / theta;
            double fact = std::max(q_norm / (theta * v_max), Fr / fr_max);
            if (fact > 1.0)
            {
                limited_cell_idx.push_back(idx_cell);
                thetaW_velocity_x[idx_cell] = q / fact;
                thetaW_velocity_y[idx_cell] = r / fact;
            }
        }
    }
    //==== Eliminate abnormally large velocities
    return EXIT_SUCCESS;
}


bool swDdp::execute(void)
/**
 * @brief Main loop for the SWES equations
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 * - all attributes (main loop)
 */
{
    ///-- Set numerical parameters locally
    t_cpu_beg = clock();
    total_t_clock_hyperb = small_real;
    total_t_clock_hyperb_lrstates = small_real;
    total_t_clock_hyperb_fluxes = small_real;
    total_t_clock_hyperb_balance = small_real;
    total_t_clock_hyperb_divcorr = small_real;
    total_t_clock_hyperb_timestep = small_real;
    total_t_clock_parab = small_real;
    total_t_clock_source = small_real;
    total_t_clock_storage = small_real;
    

    //== set numerical parameters locally
    t = t_0;
    double t_1 = clock();
    storeResultMaps(map_storage_mode, t);
    if (probes_nb > 0)
    {
        storeResultProbes(t);
    }
    total_t_clock_storage += clock() - t_1;


    ///- Time loop
    while (t < t_max)
    {
        ///-- get BC
        computeBC(t);
        //== get BC
        ///-- Hyperbolic step (always executed - should it become conditional?)
        t_1 = clock();
        if (hypCompute(div_corr_option, n_it_div_max, num_scheme, dt) != EXIT_SUCCESS)
        {
            logError() << "Error while computing hyperbolic part";
            return EXIT_FAILURE;
        }

        total_t_clock_hyperb += clock() - t_1;
        //== Hyperbolic step (always executed - should it become conditional?)

        ///-- Parabolic step (optional)
        t_1 = clock();
        if (model_momentum_diffusion != "none")
        {
            parabComputeExplicit();
        }
        total_t_clock_parab += clock() - t_1;
        //== Parabolic step (optional)

        ///-- Source term: friction (optional)
        t_1 = clock();
        /// @todo Should it be tested this way or within a source term collection ? Decide.
        if (model_friction != "none")
        {
            sourceFriction(model_friction, dt);
        }
        total_t_clock_source += clock() - t_1;
        //== Source term: friction (optional)

        ///- Source term: singular head loss (if applicable)
        if (model_singular_head_loss != "none")
        {
            sourceSingularHeadLoss(dt);
        }
        //== Source term: singular head loss (if applicable)

        ///- Runoff generation (if applicable)
        if (model_runoff != "none")
        {
            if (model_porosity_rainfall == "everywhere")
            {
                sourceRunoffEverywhere();
            }
            else if (model_porosity_rainfall == "water_only")
            {
                sourceRunoffWateronly();
            }
        }
        //== Runoff generation (if applicable)

        ///- Infiltrataion (if applicable)
        if (model_infiltration != "none")
        {
            sourceInfiltration();
        }
        //== Infiltration (if applicable)

        ///-- Source term: wind(optional)
        if (model_wind != "none")
        {
            sourceWind();
        }
        //==  Source term: wind(optional)
        if (correctVelocity() != EXIT_SUCCESS)
        {
            logError() << "Error while correcting velocity field";
            return EXIT_FAILURE;
        }

        ///-- Time incrementation
        t += dt;
        //== Time incrementation

        ///-- Store results if needed
        t_1 = clock();
        if (store_map)
        {
            std::cout << "\r  ----> calling storeResultMaps (" << t << ")" << std::flush;
            storeResultMaps(map_storage_mode, t);
        }
        if ((store_time_series) && (probes_nb > 0))
        {
            std::cout << "\r  ----> calling storeResultTimeSeries (" << t << ")" << std::flush;
            storeResultProbes(t);
        }
        total_t_clock_storage += clock() - t_1;
        //== Store results if needed
    }
    //== Time loop
    t_1 = clock();
    storeResultMaps(map_storage_mode, t);
    
    if (probes_nb > 0)
    {
        storeResultProbes(t);
        probes_output_stream.close();
    }
    total_t_clock_storage += clock() - t_1;

    double t_cpu_end = clock();

    logEmpt();
    logEmpt();
    logInfo() << "Total clock time = " << (t_cpu_end - t_cpu_beg) / CLOCKS_PER_SEC;
    logInfo() << "Hyperbolic time  = " << total_t_clock_hyperb / CLOCKS_PER_SEC << " ("
              << int(100 * total_t_clock_hyperb / (t_cpu_end - t_cpu_beg)) << " % )";
    logInfo() << "   - LRStates    = " << total_t_clock_hyperb_lrstates / CLOCKS_PER_SEC << " ("
        << int(100 * total_t_clock_hyperb_lrstates / (t_cpu_end - t_cpu_beg)) << " % )";
    logInfo() << "   - TimeStep    = " << total_t_clock_hyperb_timestep / CLOCKS_PER_SEC << " ("
        << int(100 * total_t_clock_hyperb_timestep / (t_cpu_end - t_cpu_beg)) << " % )";
    logInfo() << "   - Fluxes      = " << total_t_clock_hyperb_fluxes / CLOCKS_PER_SEC << " ("
        << int(100 * total_t_clock_hyperb_fluxes / (t_cpu_end - t_cpu_beg)) << " % )";
    logInfo() << "   - Balance     = " << total_t_clock_hyperb_balance / CLOCKS_PER_SEC << " ("
        << int(100 * total_t_clock_hyperb_balance / (t_cpu_end - t_cpu_beg)) << " % )";
    logInfo() << "   - Div. corr.  = " << total_t_clock_hyperb_divcorr / CLOCKS_PER_SEC << " ("
        << int(100 * total_t_clock_hyperb_divcorr / (t_cpu_end - t_cpu_beg)) << " % )";
    logInfo() << "Parabolic time   = " << total_t_clock_parab / CLOCKS_PER_SEC << " ("
        << int(100 * total_t_clock_parab / (t_cpu_end - t_cpu_beg)) << " % )";
    logInfo() << "Source time      = " << total_t_clock_source / CLOCKS_PER_SEC << " ("
              << int(100 * total_t_clock_source / (t_cpu_end - t_cpu_beg)) << " % )";
    logInfo() << "Storage time     = " << total_t_clock_storage / CLOCKS_PER_SEC << " ("
              << int(100 * total_t_clock_storage / (t_cpu_end - t_cpu_beg)) << " % )";
    logEmpt();
    return EXIT_SUCCESS;
}


/** @todo: adapt code */
void swDdp::importAttributesFromSolverParam(swSolverParam *solver_param)
/**
 * @brief sets class attributes by looking into m_solver_param
 * @param[in] m_solver_param: solverParam object containing all solver parameters
 * @brief modifies required (local) computational options
 * @par Modified attributes:
 * - all swSwes.h dedicated attributes
 */
{
    biblio_reference_list = solver_param->getBiblioReferenceList(); 
    cfl_max = solver_param->getCFLMax();
    chi_phi_max = solver_param->getDdpChiPhiMax();
    chi_theta_max = solver_param->getDdpChiThetaMax();
    div_corr_option = solver_param->getDivergenceCorrection();
    dt_map = solver_param->getDtMap();
    dt_max = solver_param->getDtMax();
    dt_time_series = solver_param->getDtTimeSeries();
    fr_max = solver_param->getFrMax();
    h_min = solver_param->getHMin();
    h_min_ddp = solver_param->getHMinDdp();
    log_write = solver_param->getLogWrite();
    map_storage_mode = solver_param->getStoreFlag();
    model_friction = solver_param->getModelFriction();
    model_hyp = solver_param->getModelHyp();
    model_infiltration = solver_param->getModelInfiltration();
    model_momentum_diffusion = solver_param->getModelMomentumDiffusion();
    model_momentum_dissipation = solver_param->getModelMomentumDissipation();
    model_porosity_rainfall = solver_param->getModelPorosityRainfall();
    model_runoff = solver_param->getModelRunoff();
    model_singular_head_loss = solver_param->getModelSingularHeadLoss();
    model_wind = solver_param->getModelWind();
    momentum_diff_max_iter = solver_param->getMomentumDiffMaxIter();
    momentum_diff_convergence_eps = solver_param->getMomentumDiffConvergenceEps();
    negative_depth_detection = solver_param->getNegativeDepthDetection();
    num_scheme = solver_param->getNumScheme();
    n_decimal = solver_param->getNDecimal();
    n_it_div_max = solver_param->getNItDivMax();
    n_por_tab_max = solver_param->getNPorTabMax();
    output_format = solver_param->getOutputFormat();
    probe_names = solver_param->getProbeNames();
    probe_x = solver_param->getProbeX();
    probe_y = solver_param->getProbeY();
    topo_discretization = solver_param->getTopoDiscretization();
    t_0 = solver_param->getT0();
    t_max = solver_param->getTMax();
    v_max = solver_param->getVMax();
}


bool swDdp::init(swMesh *mesh)
/**
 * @brief Initializes the SWES model
 * @param[in] mesh: contains all mesh-related fields
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 */
{
    logInfo() << "swDdp::init - START";

    ///-- Prepare output binary file for maps
    if ((output_format == "bin") || (output_format == "all"))
    {
        std::string filename = m_io->absolutePath(swFile::output_map_binary_file);
        std::string geo_filename = m_io->absolutePath(swFile::input_geo_file);
        std::vector<std::string> param_name = {};
        if (map_storage_mode == "default")
        {
            param_name = {"water_depth", "unit_discharge_x", "unit_discharge_y", "unit_disch_norm" };
        }
        else if (map_storage_mode == "all")
        {
            param_name = {"water_depth", "unit_discharge_x", "unit_discharge_y", "unit_disc_norm",
                          "surface_elevation", "velocity_x", "velocity_y", "velocity_norm",
                          "froude", "specific_head", "total_head",
                          "phi_by_wat_depth", "phi_by_disc_x", "phi_by_disc_y", "phi_by_disc_norm"};
        }
        if (swIO::initBinFileMapDecimal(filename, geo_filename, param_name, t_0, t_max, dt_map) != EXIT_SUCCESS)
        {
            logError() << "Error while initializing the binary output file";
            return EXIT_FAILURE;
        }
    }
    //==  Prepare output binary file for maps

    ///-- Prepare output file for probes (first lines contain probes name and location)
    probes_nb = probe_names.size();
    probes_filename = m_io->absolutePath(swFile::output_probes_file);
    if (probes_nb > 0)
    {
        probes_output_stream.open(probes_filename);
        for (size_t p = 0 ; p < probes_nb ; p++)
        {
            std::string name = probe_names[p];
            double x = probe_x[p];
            double y = probe_y[p];
            swCell* probe_cell = mesh->getCellFromCoordinates(x,y);
            if (probe_cell == nullptr)
            {
                logError() << "Error while looking for probe coordinates. Probe could be out of the mesh";
                return EXIT_FAILURE;
            }
            probe_cell_idx.push_back(probe_cell->getIdx());
            std::string probes_file_line = name + "\t" + swString::numberToStr(x, 'e', 10) + "\t" + swString::numberToStr(y, 'e', 10) + "\n";
            swIO::writeTxtLine(probes_output_stream, probes_file_line);
        }
        std::string probes_file_line = "time\t";
        for (size_t p = 0 ; p < probes_nb ; p++)
        {
            std::string name = probe_names[p];
            probes_file_line += name + "_zb\t" + name + "_h\t" + name + "_q\t" + name + "_r\t";
        }
        probes_file_line += "\n";
        swIO::writeTxtLine(probes_output_stream, probes_file_line);
    }
    //==  Prepare output file for probes (first lines contain probes name and location)
    

    ///-- Initialize mesh   
	logInfo() << "swDdp::init - get mesh information() " ;
    setMeshTables(mesh);
	//== Initialize mesh

    ///-- Initialize variables & parameters 

    ///------Porosity (compulsory)
    setPorosityParameters();
    //=======Porosity (compulsory)

	///---- Set up initial conditions 
    //-- NB: variables must be declared in .h
    setInitialConditions();
    //==== Set up initial conditions


    ///---- Set up boundary conditions
    logInfo() << "Setting up boundary conditions";
    if (swParent::initBoundaryConditions() != EXIT_SUCCESS)
    {
        logError() << "Error while initializing boundary conditions";
        return EXIT_FAILURE;
    }
    logInfo() << "  ...done";
    //==== Set up boundary conditions


    ///---- Set up parabolic model (if applicable)
    if (model_momentum_diffusion != "none" || model_transport_diffusion != "none")
    {
        setParabParameters();
    }
    ///---- Set up parabolic (if applicable)

    ///----Read parameter files

    ///------Boussinesq (compulsory)
    setBoussinesqParameter();
    //=======Boussinesq (compulsory)
    
    ///------Friction
    logInfo() << "swDdp::init - friction_model = " << model_friction;
    if (model_friction != "none")
    {
        setFrictionParameters();
    }
    //====== Friction
        
    ///------ Momentum dissipation
    logInfo() << "swDdp::init - model_momentum_dissipation = " << model_momentum_dissipation;
    if (model_momentum_dissipation != "none")
    {
        setMomentumDissipationParameters();
    }
    //====== Momentum dissipation

    ///------Singular head loss
    if (model_singular_head_loss != "none")
    {
        setSingularHeadLossParameters();
    }
    //=======Singular head loss

    ///------Runoff
    if (model_runoff != "none")
    {
        setRunoffParameters();
        setPrecipitationStationCodes();
        setPrecipitationTimeSeriesTable();
    }
    //=======Runoff

    ///------Infiltration
    if (model_infiltration != "none")
    {
        setInfiltrationParameters();
    }
    //=======Infiltration

    ///------Wind
    if (model_wind != "none")
    {
        setWindStationCodes();
        setWindVelocityTimeSeriesTables();
    }
    //=======Wind

    ///---- Set up transport model (if applicable)
    if (model_transport_diffusion != "none")
    {
        setTransportDiffusionParameters();
    }


    //==== Read parameter files
    //== Initialize variables & parameters

	///-- Compute edge factors 
    logInfo() << "swDdp::init - compute edge factors";
    computeEdgeFactors();
    //== Compute edge factors 

    ///-- Resize fluxes and Riemann states
    lambda_max_bnd.resize(edge_bnd_nb);
    lambda_max_int.resize(edge_int_nb);
    lambda_minus_bnd.resize(edge_bnd_nb);
    lambda_minus_int.resize(edge_int_nb);
    lambda_plus_bnd.resize(edge_bnd_nb);
    lambda_plus_int.resize(edge_int_nb);
    theta_riemann_bnd.resize(edge_bnd_nb);
    theta_riemann_int.resize(edge_int_nb);
    theta_u_riemann_bnd.resize(edge_bnd_nb);
    theta_u_riemann_int.resize(edge_int_nb);
    theta_v_riemann_bnd.resize(edge_bnd_nb);
    theta_v_riemann_int.resize(edge_int_nb);
    z_riemann_bnd.resize(edge_bnd_nb);
    z_riemann_int.resize(edge_int_nb);
    F1_bnd.resize(edge_bnd_nb);
    F2_bnd.resize(edge_bnd_nb);
    F3_bnd.resize(edge_bnd_nb);
    F1_int.resize(edge_int_nb);
    F2_int.resize(edge_int_nb);
    F3_int.resize(edge_int_nb);
    s2L_int.resize(edge_int_nb);
    s2R_bnd.resize(edge_bnd_nb);
    s2R_int.resize(edge_int_nb);
    //== Resize fluxes and Riemann states


    //-- Constistency check
    checkParameterConsistency();
    //== Constistency check

    //-- Write Log files (if applicable)
    if (log_write)
    {
        logInfo() << "swDdp::init - write log files";
        writeLogs();
    }
    //== Write Log files (if applicable)

    logInfo() << "swDdp::init - done \n";
    return EXIT_SUCCESS;
}


double swDdp::lambdaMinus(double const& u_cell, double const& c, double const& x_phi, double const& x_theta, double const& b)
/**
 * @brief Computes lambda Minus
 * @param u_cell: velocity in the cell
 * @param c: wave speed sqrt(g.h)
 * @param x_phi: phiW / phiG
 * @param x_theta: thetaW / thetaG
 * @param b: Boussinesq coefficient
 * @return lambda_minus
 */
{
    return b * x_theta * u_cell - std::sqrt(x_theta * x_theta * (b - 1./std::max(x_phi,small_real)) * b * u_cell * u_cell + c * c / std::max(x_theta, small_real));
}


double swDdp::lambdaPlus(double const& u_cell, double const& c, double const& x_phi, double const& x_theta, double const& b)
/**
 * @brief Computes lambda Plus
 * @param u_cell: velocity in the cell
 * @param c: wave speed sqrt(g.h)
 * @param x_phi: phiW / phiG
 * @param x_theta: thetaW / thetaG
 * @param b: Boussinesq coefficient
 * @return lambda_plus
 */
{
    return b * x_theta * u_cell + std::sqrt(x_theta * x_theta * (b - 1./std::max(x_phi,small_real)) * b * u_cell * u_cell + c * c / std::max(x_theta, small_real));
}


double swDdp::phiGBndOfZ(const size_t& idx_edge_bnd, const double& z)
/**
 * @brief Compute phi(z) in any given boundary edge, where phi is tabulated from the porosity data
 * @param idx_edge_bnd: boundary edge index in which the porosity is computed
 * @param z: altitude where the boundary edge porosity is evaluated
 * @return evaluated porosity in boundary edge idx_edge_bnd, at altitude z
 */
{
    std::vector<double> z_tab = edge_bnd_z_tab[idx_edge_bnd];
    if (z <= z_tab[0])
    {
        return 0.;
    }
    else
    {
        size_t level = 0;
        while ((level < z_tab.size()) && (z > z_tab[level]))
        {
            ++level;
        }
        return phiG_bnd_tab[idx_edge_bnd][level - 1];
    }
}


double swDdp::phiGIntOfZ(const size_t& idx_edge_int, const double& z)
/**
 * @brief Compute phi(z) in any given internal edge, where phi is tabulated from the porosity data
 * @param idx_edge_int: internal edge index in which the porosity is computed
 * @param z: altitude where the internal edge porosity is evaluated
 * @return evaluated porosity in internal edge idx_edge_int, at altitude z
 */
{
    std::vector<double> z_tab = edge_int_z_tab[idx_edge_int];
    if (z <= z_tab[0])
    {
        return 0.;
    }
    else
    {
        size_t level = 0;
        while ((level < z_tab.size()) && (z > z_tab[level]))
        {
            ++level;
        }
        return phiG_int_tab[idx_edge_int][level - 1];
    }
}


double swDdp::phiWOfZ(const size_t& idx_cell, const double& z)
/**
 * @brief Compute phi(z) in any given cell, where phi is tabulated from the porosity data
 * @param idx_cell: cell index in which the porosity is computed
 * @param z: altitude where the cell porosity is evaluated
 * @return evaluated porosity in cell idx_cell, at altitude z
 */
{
    std::vector<double> z_tab = cell_z_tab[idx_cell];
    if (z <= z_tab[0])
    {
        return 0.;
    }
    else
    {
        size_t level = 0;
        while ((level < z_tab.size()) && (z > z_tab[level]))
        {
            ++level;
        }
        return phiW_tab[idx_cell][level - 1];
    }
}


size_t swDdp::porosityLineSize(const size_t& line_nb, const double& porosity_law_type, const std::vector<std::string> line_txt)
/**
 * @brief Compute expected number of parameters in each line of DDP porosity files as a function of the porosity law type. 
 * The returned number includes the slot for the first value, i.e. porosity_law_type.
 * @param line_nb: current line
 * @param porosity_law_type: type of porosity law, defining the number of parameters expected
 * @param line_txt: content of the whole line, size should be compatible with the value of porosity_law_type
 * @returns number of parameters that should be read
 */
{
    size_t line_size = line_txt.size();
    if ((porosity_law_type == -2) || (porosity_law_type == -1)) // Arbitrary number n of couples (phi, z)
    // Line should be: law_type | n | z_0 | phi_0 | z_1 | phi_1 | ... | z_{n-1} | phi_{n-1}
    // There are 2n+2 records in the line
    {
        bool conversion_success = true;
        int n_couples = swString::stringToInt(line_txt[1], conversion_success);
        if (conversion_success == false)
            {
                logError() << "line is misformatted: format conversion";
                logFatal(__FILE__, __func__, __LINE__);
            }
        if (line_size == 2*n_couples + 2)
        {
            return line_size;
        }
        else
        {
            logError() << "line " << line_nb << " is misformatted";
            logError() << "  porosity_law_type = " << porosity_law_type << " is not compatible with " << line_size-1 << "parameters.";
            logFatal(__FILE__, __func__, __LINE__);
            return 0;
        }
    }
    else if (porosity_law_type == 0) // Arbitrary number n of (increasing) z levels, porosity to be interpolated between z_0 and z_{n-1}
    // Line should be: 0 | n | z_0 | z_1 | ... | z_{n-1}
    // There are n+2 records in the line
    {
        bool conversion_success = true;
        int n_levels = swString::stringToInt(line_txt[1], conversion_success);
        if (conversion_success == false)
            {
                logError() << "line is misformatted: format conversion";
                logFatal(__FILE__, __func__, __LINE__);
            }
        if (line_size == n_levels + 2)
        {
            return line_size;
        }
        else
        {
            logError() << "line " << line_nb << " is misformatted";
            logError() << "  porosity_law_type = " << porosity_law_type << " is not compatible with " << line_size-1 << "parameters.";
            logFatal(__FILE__, __func__, __LINE__);
            return 0;
        }
    }
    else if (porosity_law_type == 1) // Uniform case, 2 parameters expected: z_0 and phi
    // Line should be: 1 | 2 | z_0 | phi
    // There are 4 records in the line
    {
        if  (line_size == 4)
        {
            return line_size;
        }
        else
        {
            logError() << "line " << line_nb << " is misformatted";
            logError() << "  porosity_law_type = 1, uniform case expecting 3 more parameters (2, z_0 and phi), but got " << line_size-1 << ".";
            logFatal(__FILE__, __func__, __LINE__);
            return 0;
        }
    }
    else if (porosity_law_type == 2) // Linear case, 4 parameters expected: z0, z1, phi0, phi1
    // Line should be: 2 | 4 | z_0 | phi_0 | z_1 | phi_1
    // There are 6 records in the line
    {
        if  (line_size == 6)
        {
            return line_size;
        }
        else
        {
            logError() << "line " << line_nb << " is misformatted";
            logError() << "  porosity_law_type = 2, linear case expecting 5 more parameters (4, z0, z1, phi0, phi1), but got " << line_size-1 << ".";
            logFatal(__FILE__, __func__, __LINE__);
            return 0;
        }
    }
    else if (porosity_law_type == 3) // Linear case + discontinuity value, 5 parameters expected: z0, z1, phi0, phi1, phi2
    // Line should be: 3 | 5 | z_0 | phi_0 | z_1 | phi_1 | phi_2
    // There are 7 records in the line
    {
        if  (line_size == 7)
        {
            return line_size;
        }
        else
        {
            logError() << "line " << line_nb << " is misformatted";
            logError() << "  porosity_law_type = 3, case expecting 5 more parameters (5, z0, z1, phi0, phi1, phi2), but got " << line_size-1 << ".";
            logFatal(__FILE__, __func__, __LINE__);
            return 0;
        }
    }
    else 
    {
        logError() << "line " << line_nb << " is misformatted";
        logError() << "  porosity_law_type = " << porosity_law_type << " is not compatible with " << line_size-1 << "parameters.";
        logFatal(__FILE__, __func__, __LINE__);
        return 0;
    }
}


double swDdp::pressureForceGammaOfZ(const size_t& idx_edge_int, const double& z)
/**
 * @brief Compute the pressure force against PhiG from z
 * @param idx_edge_int: internal edge index
 * @param z: free surface elevation
 * @return
 */
{
    double pressure_force = 0;
    size_t z_tab_size = edge_int_z_tab[idx_edge_int].size();
    if (z > edge_int_z_tab[idx_edge_int][0])
    {
        size_t level = 0;
        while ((level < z_tab_size - 1) &&
               (z >= edge_int_z_tab[idx_edge_int][level+1]))
        {
            double z1 = edge_int_z_tab[idx_edge_int][level];
            double z2 = edge_int_z_tab[idx_edge_int][level+1];
            double dz1 = z - (z1 + z2) / 2.;
            double dz2 = z2 - z1;
            pressure_force += phiG_int_tab[idx_edge_int][level] * dz1 * dz2;
            ++level;
        }
        double dz1 = std::max(0., z - edge_int_z_tab[idx_edge_int][level]);
        pressure_force += phiG_int_tab[idx_edge_int][level] * dz1 * dz1 / 2.;
    }
    pressure_force *= gravity;
    return pressure_force;
}


double swDdp::pressureForceOmegaOfZ(const size_t& idx_cell, const double& z)
/**
 * @brief Compute the pressure force against PhiW from z
 * @param idx_cell: cell index
 * @param z: free surface elevation
 * @return
 */
{
    double pressure_force = 0;
    size_t z_tab_size = cell_z_tab[idx_cell].size();
    if (z > cell_z_tab[idx_cell][0])
    {
        size_t level = 0;
        while ((level < z_tab_size - 1) &&
               (z >= cell_z_tab[idx_cell][level+1]))
        {
            double z1 = cell_z_tab[idx_cell][level];
            double z2 = cell_z_tab[idx_cell][level+1];
            double dz1 = z - (z1 + z2) / 2.;
            double dz2 = z2 - z1;
            pressure_force += phiW_tab[idx_cell][level] * dz1 * dz2;
            ++level;
        }
        double dz1 = std::max(0., z - cell_z_tab[idx_cell][level]);
        pressure_force += phiW_tab[idx_cell][level] * dz1 * dz1 / 2.;
    }
    pressure_force *= gravity;
    return pressure_force;
}


double swDdp::pressureWaveCelerity(double const &theta, double const &phi)
/**
 * @brief 
 * 
 * @param theta 
 * @param phi 
 * @return double 
 */
{
    return sqrt(gravity * theta / std::max(phi, small_real));
}


bool swDdp::readPorosityFile(const std::string& file_name, const size_t& record_nb, std::vector<std::string>& headers, std::vector<std::vector<double>>& values)
/**
 * @brief Specific routine to read porosity file (dedicated format)
 * @param[in] file_name : name of the file to be read
 * @param[in] record_nb : number of records (number of lines in the file if parameter is not uniform, could be cells or edges)
 * @param[out] headers [n] : name of the n-th parameter in the file
 * @param[out] values [k][n] : the numerical value of the k-th record of the n-th parameter in the file
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 */
{
    //---- read DDP porosity file 
    std::ifstream file_r;
    file_r.open(file_name);
    if (file_r.is_open() == false)
    {
        logError() << "    Unable to open file " << file_name;
        return EXIT_FAILURE;
    }

    logInfo() << "Importing porosity parameters from file:" << file_name;

    size_t line_nb = 0;
    std::string line_txt;

    //------ Read distribution type (uniform or not)
    line_txt = swIO::readLine(file_r);
    ++line_nb;
    std::vector<std::string> line_split = swIO::stringSplit(line_txt);
    if (line_split.size() != 2)
    {
        logError() << "line " << line_nb << " is misformatted: number of arguments = " << line_split.size();
        return EXIT_FAILURE;
    }

    bool conversion_success = true;
    const int unif = swString::stringToInt(line_split[1], conversion_success);
    if (conversion_success == false)
    {
        logError() << "line is misformatted: format conversion";
        swIO::errorReadingLine(file_r, file_name, line_txt, line_nb, __LINE__, __func__, __FILE__);
    }

    if (!((unif == 0) || (unif == 1)))
    {
        logError() << "Incorrect value for 'Unif' parameter in porosity file";
        logError() << "  allowed value for parameter 'Unif' : 0 or 1";
        return EXIT_FAILURE;
    }
    logInfo() << "Porosity parameter uniformity = " << unif;
    //====== Read distribution type (uniform or not)


    //------ Read dummy line ("==== Default Param")
    line_txt = swIO::readLine(file_r);
    ++line_nb;
    //====== Read dummy line ("==== Default Param")

    //------ Read header (name of variables)
    line_txt = swIO::readLine(file_r);
    ++line_nb;
    line_split = swIO::stringSplit(line_txt);
    
    /// The header line format is not known, since the number of parameters may change from one record (cell) to the other
    //-- Prepare array
    headers.resize(line_split.size());
    values.resize(record_nb);
    //== Prepare arrays
    
    //-- initialize vector for storing header
    for (size_t i = 0; i < line_split.size(); i++)
    {
        headers[i] = line_split[i];
        logInfo() << "  param number and name = " << i << ", " << headers[i];
    }
    //== initialize vector for storing header
    //======= Read header (name of variables)

    //------ Uniform case
    if (unif == 1)
    {
        line_txt = swIO::readLine(file_r);
        ++line_nb;
        line_split = swIO::stringSplit(line_txt);
        

        //-- Manage porosity law cases
        double porosity_law_type = swString::stringToInt(line_split[0], conversion_success);
        if (conversion_success == false)
        {
            logError() << "line is misformatted: format conversion";
            swIO::errorReadingLine(file_r, file_name, line_txt, line_nb, __LINE__, __func__, __FILE__);
        }
        size_t n_params = porosityLineSize(line_nb, porosity_law_type, swIO::stringSplit(line_txt));
        //== Manage porosity law cases

        //-- Fill porosity array
        //---- First record
        values[0] = {};
        for (size_t idx_param = 0; idx_param < n_params; ++idx_param) 
        {
            bool conversion_success = true;
            double val = swString::stringToDouble(line_split[idx_param], conversion_success);
            logInfo() << "Assigning uniform variable - record, value = " << idx_param << val;
            if (conversion_success == false)
            {
                logError() << "line is misformatted: format conversion";
                swIO::errorReadingLine(file_r, file_name, line_txt, line_nb, __LINE__, __func__, __FILE__);
                return EXIT_FAILURE;
            }
            else
            {
                values[0].push_back(val);
            }
        }
        //==== First record
        //---- Other records 
        for (size_t idx_record = 1; idx_record < record_nb; idx_record++)
        {
            values[idx_record] = {};
            for (size_t idx_param = 0; idx_param < n_params; ++idx_param)
            {
                values[idx_record].push_back(values[0][idx_param]);
            }
        }
        //==== Other records
        //== Fill porosity array
        file_r.close();
    }
    //====== Uniform case
    // 
    // 
    //------ Variable case
    else
    {
        /* Read 2 commented lines */
        line_txt = swIO::readLine(file_r);// Default values (not used here)
        ++line_nb;
        line_txt = swIO::readLine(file_r);// ==== Distrib
        ++line_nb;

        const size_t nb_dummy_lines = 5;
        for (size_t idx_record = 0; idx_record < record_nb; idx_record++)
        {
            if (file_r.eof() == true)
            {
                file_r.close();
                logError() << "Not enough values in file" << file_name;
                logError() << "Expected/existing number of lines containing porosity parameters:"
                    << record_nb << line_nb - nb_dummy_lines - 1;
                return EXIT_FAILURE;
            }

            line_txt = swIO::readLine(file_r);
            ++line_nb;
            line_split = swIO::stringSplit(line_txt);

            //-- Manage porosity law cases
            double porosity_law_type = swString::stringToInt(line_split[0], conversion_success);
            size_t n_params = porosityLineSize(line_nb, porosity_law_type, swIO::stringSplit(line_txt));
            //== Manage porosity law cases
            
            //-- Fill porosity array
            values[idx_record] = {};
            for (size_t idx_param = 0; idx_param < n_params; idx_param++)
            {
                bool conversion_success = true;
                double val = swString::stringToDouble(line_split[idx_param], conversion_success);
                if (conversion_success == false)
                {
                    logError() << "line is misformatted: format conversion";
                    swIO::errorReadingLine(file_r, file_name, line_txt, line_nb, __LINE__, __func__, __FILE__);
                    return EXIT_FAILURE;
                }
                values[idx_record].push_back(val);
            }
            //== Fill porosity array
            ++line_nb;
        }
        //====== Variable case

        // Warning if there are additional records in file
        if (!file_r.eof())
        {
            logWarn() << "Additional lines found in file " << file_name << " will be ignored.";
        }
        // Close file
        file_r.close();
    }
    logInfo() << "...done.";
    //==== read DDP porosity file file
    return EXIT_SUCCESS;
}


/** @todo: adapt code */
void swDdp::setBoussinesqParameter(void)
/**
 * @brief Reads the distribution of the Boussinesq coefficient from file
 * @param Nothing
 * @return Nothing
 * @par Modified attributes:
 *  - boussinesq
 */
{
    logInfo() << "swDdp::init - reading Boussinesq coeff. map";
    //-- Set up pre-defined input file name as define in file.h
    std::string file_name;
    file_name = swFile::input_boussinesq_file;
    file_name = m_io->absolutePath(file_name);
    std::vector<std::string> variable_name;
    std::vector<std::vector<double>> param_array;
    swIO::readVectorParameterFile(file_name, 1, cell_nb, variable_name, param_array);
    boussinesq.resize(cell_nb);
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        boussinesq[idx_cell] = param_array[0][idx_cell];
    }
    logInfo() << "           - Boussinesq coeff. map read \n";
}


/** @todo: adapt code */
void swDdp::setFrictionParameters(void)
/**
 * @brief read information from friction file (could be Chez, Manning, Strickler)
 * @par modified attributes:
 *  - friction_coeff_direction_1
 *  - friction_coeff_direction_2
 *  - friction_alpha
 */
{
    logInfo() << "swDdp::init - reading friction map";
    //-- Set up pre-defined input file name as define in file.h
    std::string file_name;
    if (model_friction == "chezy")
    {
        file_name = swFile::input_friction_chezy_map;
    }
    else if (model_friction == "manning")
    {
        file_name = swFile::input_friction_manning_map;
    }
    else if (model_friction == "strickler")
    {
        file_name = swFile::input_friction_strickler_map;
    }
    else
    {
        logError() << "Friction model not handled = " << model_friction;
        logFatal(__FILE__, __func__, __LINE__);
    }
    file_name = m_io->absolutePath(file_name);
    std::vector<std::string> variable_name;
    std::vector<std::vector<double>> param_array;
    swIO::readVectorParameterFile(file_name, 3, cell_nb, variable_name, param_array);
    friction_coeff_direction_1.resize(cell_nb);
    friction_coeff_direction_2.resize(cell_nb);
    friction_alpha.resize(cell_nb);
    friction_cos_alpha.resize(cell_nb);
    friction_sin_alpha.resize(cell_nb);
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        friction_coeff_direction_1[idx_cell] = param_array[0][idx_cell];
        friction_coeff_direction_2[idx_cell] = param_array[1][idx_cell];
        double alpha = param_array[2][idx_cell] * atan(1.0) / 45.0;  // Convert degrees to radians
        friction_alpha[idx_cell]     = alpha;
        friction_cos_alpha[idx_cell] = cos(alpha);
        friction_sin_alpha[idx_cell] = sin(alpha);
    }
    logInfo() << "             - friction map read \n";
}


/** @todo: adapt code */
void swDdp::setInfiltrationParameters(void)
/**
 * @brief read information from infiltration file (could be h_proportional, Horton)
 * @par modified attributes:
 *  - infiltration_param
 */
{
    logInfo() << "\n \n";
    logInfo() << "swDdp::init - reading infiltration parameter map";
    logInfo() << "             - model_infiltration = " << model_infiltration;
    std::string file_name;
    std::vector<std::string> variable_name;
    std::vector<std::vector<double>> param_array;

    infiltration_param.resize(cell_nb);
    if (model_infiltration == "h_proportional")
    {
        file_name = m_io->absolutePath(swFile::input_infiltration_h_proportional_parameter_map);
        infiltration_param.resize(cell_nb);
        swIO::readVectorParameterFile(file_name, 1, cell_nb, variable_name, param_array);
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            infiltration_param[idx_cell] = param_array[0][idx_cell];
        }
    }
    else if (model_infiltration == "horton")
    {
       file_name = m_io->absolutePath(swFile::input_infiltration_horton_parameter_map);
       infiltration_param.resize(cell_nb);
       swIO::readVectorParameterFile(file_name, 1, cell_nb, variable_name, param_array);
       for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
       {
           infiltration_param[idx_cell] = param_array[0][idx_cell];
       }
    }
    else
    {
        logError() << "Infiltration model not handled = " << model_infiltration;
        logFatal(__FILE__, __func__, __LINE__);
    }
    logInfo() << "             - Infiltration parameter map read \n";
}


/** @todo: adapt code */
void swDdp::setInitialConditions(void)
/**
 * @brief read information from initial conditions file and assign it to state variables
 * @par modified attributes:
 *  - water_depth
 *  - unit_discharge_x
 *  - unit_discharge_y
 */
{
    logInfo() << "swDdp::init - reading initial conditions";
    std::string file_name = m_io->absolutePath(swFile::input_initial_conditions_file);

    std::vector<std::string> variable_name;
    std::vector<std::vector<double>> param_array;
    swIO::readVectorParameterFile(file_name, 3, cell_nb, variable_name, param_array);

    thetaW.resize(cell_nb);
    thetaW_velocity_x.resize(cell_nb);
    thetaW_velocity_y.resize(cell_nb);
    thetaW_old.resize(cell_nb);
    thetaW_velocity_x_old.resize(cell_nb);
    thetaW_velocity_y_old.resize(cell_nb);

    ///------ Assign initial water depth (if h or z present in variable names)
    bool found = false;
    auto it = std::find(variable_name.begin(), variable_name.end(), "h");
    size_t column = it - variable_name.begin();
    if (column < 3) /* Index has been found before end of vector*/
    {
        found = true;
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            double zs = param_array[column][idx_cell] + cell_z_tab[idx_cell][0];
            thetaW[idx_cell] = thetaWOfZ(idx_cell, zs);
        }
    }
    it = std::find(variable_name.begin(), variable_name.end(), "z");
    column = it - variable_name.begin();
    if (column < 3) /* Index has been found before end of vector*/
    {
        found = true;
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            double zs = param_array[column][idx_cell];
            thetaW[idx_cell] = thetaWOfZ(idx_cell, zs);
        }
    }
    if (!found)
    {
        logError() << "Could find neither h nor z in initial conditions. Stopping here.";
        logFatal(__FILE__, __func__, __LINE__);
    }
    //====== Assign initial water depth (if h or z present in variable names)

    ///------ Assign initial x-unit discharge (if q or u present in variable names)
    found = false;
    it = std::find(variable_name.begin(), variable_name.end(), "q");
    column = it - variable_name.begin();
    if (column < 3) /* Index has been found before end of vector*/
    {
        found = true;
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            double q = param_array[column][idx_cell];
            double h = std::max(small_real, zOfThetaW(idx_cell, thetaW[idx_cell]) - cell_z_tab[idx_cell][0]);
            thetaW_velocity_x[idx_cell] = thetaW[idx_cell] * q / h;
        }
    }
    it = std::find(variable_name.begin(), variable_name.end(), "u");
    column = it - variable_name.begin();
    if (column < 3) /* Index has been found before end of vector*/
    {
        found = true;
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            double u = param_array[column][idx_cell];
            thetaW_velocity_x[idx_cell] = thetaW[idx_cell] * u;
        }
    }
    if (!found)
    {
        logError() << "Could find neither q nor u in initial conditions. Stopping here.";
        logFatal(__FILE__, __func__, __LINE__);
    }
    //====== Assign initial x-unit discharge (if q or u present in variable names)

    ///------ Assign initial y-unit discharge (if r or v present in variable names)
    found = false;
    it = std::find(variable_name.begin(), variable_name.end(), "r");
    column = it - variable_name.begin();
    if (column < 3) /* Index has been found before end of vector*/
    {
        found = true;
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            double r = param_array[column][idx_cell];
            double h = std::max(small_real, zOfThetaW(idx_cell, thetaW[idx_cell]) - cell_z_tab[idx_cell][0]);
            thetaW_velocity_y[idx_cell] = thetaW[idx_cell] * r / h;
        }
    }
    it = std::find(variable_name.begin(), variable_name.end(), "v");
    column = it - variable_name.begin();
    if (column < 3) /* Index has been found before end of vector*/
    {
        found = true;
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            double v = param_array[column][idx_cell];
            thetaW_velocity_y[idx_cell] = thetaW[idx_cell] * v;
        }
    }
    if (!found)
    {
        logError() << "Could find neither r nor v in initial conditions. Stopping here.";
        logFatal(__FILE__, __func__, __LINE__);
    }
    //====== Assign initial y-unit discharge (if r or v present in variable names)
}


/** @todo: adapt code */
void swDdp::setMomentumDissipationParameters(void)
/**
 * @brief read information from momentum dissipation file (could be Chez, Manning, Strickler)
 * @par modified attributes:
 *  - friction_coeff_direction_1
 *  - friction_coeff_direction_2
 *  - friction_alpha
 */
{
    logInfo() << "swDdp::init - reading friction map";
    //-- Set up pre-defined input file name as define in file.h
    std::string file_name;
    file_name = swFile::input_momentum_dissipation_parameters_map;
    file_name = m_io->absolutePath(file_name);
    std::vector<std::string> variable_name;
    std::vector<std::vector<double>> param_array;
    swIO::readVectorParameterFile(file_name, 3, cell_nb, variable_name, param_array);
    cell_mu_1.resize(cell_nb);
    cell_mu_2.resize(cell_nb);
    cell_mu_alpha.resize(cell_nb);
    cell_mu_cos_alpha.resize(cell_nb);
    cell_mu_sin_alpha.resize(cell_nb);
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        cell_mu_1[idx_cell] = param_array[0][idx_cell];
        cell_mu_2[idx_cell] = param_array[1][idx_cell];
        double alpha = param_array[2][idx_cell] * atan(1.0) / 45.0;  // Convert degrees to radians
        cell_mu_alpha[idx_cell]     = alpha;
        cell_mu_cos_alpha[idx_cell] = cos(alpha);
        cell_mu_sin_alpha[idx_cell] = sin(alpha);
    }
    logInfo() << "             - friction map read \n";
}


/** @todo: adapt code */
void swDdp::setParabParameters(void)
/**
 * @brief read information from parabolic part
 * @par modified attributes:
 *  - parab_edge_int_flux_factors 
 */
{
    logInfo() << "setParabParameters - top";
    //-- Prepare coefficients for matrix setup
    for (size_t idx_edge_int = 0; idx_edge_int < edge_int_nb; idx_edge_int++)
    {
        size_t idx_celL = edge_int_celL_idx[idx_edge_int];
        size_t idx_celR = edge_int_celR_idx[idx_edge_int];
        double dx = cell_centre_x[idx_celR] - cell_centre_x[idx_celL];
        double dy = cell_centre_y[idx_celR] - cell_centre_y[idx_celL];
        double d_norm = std::abs(dx * edge_int_norm_x[idx_edge_int] + dy * edge_int_norm_y[idx_edge_int]);
        parab_edge_int_flux_factors.push_back(edge_int_length[idx_edge_int] / std::max(d_norm, small_real));
    }
    //== Prepare coefficients for matrix setup
}


bool swDdp::setPorosityCellMap(const std::vector<std::vector<double>>& records)
/**
 * @todo change into boo with EXIT_SUCCESS / FAILURE ?
 * @brief set porosity arrays (z, phi and theta) in cells
 * @param[in] records : values of law_type, z and phi to be tabulated
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par modified attributes:
 *  - cell_z_tab
 *  - phiW_tab
 *  - thetaW_tab
 */
{
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        //-- Compute cell_z_tab & phiW_tab
        double law_type = records[idx_cell][0];
        if (law_type == -2) // Arbitrary number n of couples (phi, z), these values are not re-tabulated
        // Line should be: law_type | n | z_0 | phi_0 | z_1 | phi_1 | ... | z_{n-1} | phi_{n-1}
        // There are 2n+2 records in the line
        {
            size_t n_couples = int(records[idx_cell][1]);
            logWarn() << "Using pre-tabulated porosities. Number of vertical porosity levels changed to" << n_couples;
            for (size_t z_level = 0; z_level < n_couples; z_level++)
            {
                double z = cell_z_b[idx_cell] + records[idx_cell][2 + 2*z_level];
                double phi = records[idx_cell][3 + 2*z_level];
                cell_z_tab[idx_cell].push_back(z);
                phiW_tab[idx_cell].push_back(phi);
            }
        }
        else if (law_type == -1) // Arbitrary number n of couples (phi, z), these values are re-tabulated with n_por_tab_max levels
        // Line should be: law_type | n | z_0 | phi_0 | z_1 | phi_1 | ... | z_{n-1} | phi_{n-1}
        // There are 2n+2 records in the line
        {
            std::vector<double> phi_read = {};
            std::vector<double> z_read = {};
            size_t n_levels = int(records[idx_cell][1]);
            for (size_t z_level = 0; z_level < n_levels; z_level++)
            {
                z_read.push_back(cell_z_b[idx_cell] + records[idx_cell][2 + 2*z_level]);
                phi_read.push_back(records[idx_cell][3 + 2*z_level]);
            }
            if (tabulatePorosity(n_por_tab_max, z_read, phi_read, cell_z_tab[idx_cell], phiW_tab[idx_cell]) != EXIT_SUCCESS)
            {
                logError() << "Error while tabulating porosity values with law type -1";
                return EXIT_FAILURE;
            }
        }
        else if (law_type == 0) // Arbitrary number n of (increasing) z levels, porosity to be interpolated between z_0 and z_{n-1}
        // Line should be: 0 | n | z_0 | z_1 | ... | z_{n-1}
        // There are n+2 records in the line
        {
            std::vector<double> phi_read = {};
            std::vector<double> z_read = {};
            size_t n_levels = int(records[idx_cell][1]);
            for (size_t z_level = 0; z_level < n_levels; z_level++)
            {
                z_read.push_back(cell_z_b[idx_cell] + records[idx_cell][2 + z_level]);
                phi_read.push_back(z_level * 1.0 / (n_levels-1));
            }
            if (tabulatePorosity(n_por_tab_max, z_read, phi_read, cell_z_tab[idx_cell], phiW_tab[idx_cell]) != EXIT_SUCCESS)
            {
                logError() << "Error while tabulating porosity values with law type 0";
                return EXIT_FAILURE;
            }
        }
        else if (law_type == 1) // Uniform case, 2 parameters expected: z_0 and phi
        // Line should be: 1 | 2 | z_0 | phi
        // There are 4 records in the line
        {
            double z = cell_z_b[idx_cell] + records[idx_cell][2];
            double phi = records[idx_cell][3];
            std::vector<double> phi_read = {};
            std::vector<double> z_read = {};
            z_read.push_back(z);
            z_read.push_back(z + 1.0);
            phi_read.push_back(phi);
            phi_read.push_back(phi);
            if (tabulatePorosity(n_por_tab_max, z_read, phi_read, cell_z_tab[idx_cell], phiW_tab[idx_cell]) != EXIT_SUCCESS)
            {
                logError() << "Error while tabulating porosity values with law type 1";
                return EXIT_FAILURE;
            }
            phiW_tab[idx_cell][0] = phi;
        }
        else if (law_type == 2) // Linear case, 4 parameters expected: z0, z1, phi0, phi1
        // Line should be: 2 | 4 | z_0 | phi_0 | z_1 | phi_1
        // There are 6 records in the line
        {
            double z0 = cell_z_b[idx_cell] + records[idx_cell][2];
            double z1 = cell_z_b[idx_cell] + records[idx_cell][4];
            double phi0 = records[idx_cell][3];
            double phi1 = records[idx_cell][5];
            double dz = (z1 - z0) / (n_por_tab_max - 1);
            double dphi = (phi1 - phi0) / (n_por_tab_max - 1);
            cell_z_tab[idx_cell].resize(n_por_tab_max);
            phiW_tab[idx_cell].resize(n_por_tab_max);
            
            for (size_t z_level = 0; z_level < n_por_tab_max - 1; ++z_level)
            {
                cell_z_tab[idx_cell][z_level] = z0 + z_level * dz;
                phiW_tab[idx_cell][z_level] = phi0 + (z_level + 0.5) * dphi;
            }
            cell_z_tab[idx_cell][n_por_tab_max - 1] = z1;
            phiW_tab[idx_cell][n_por_tab_max - 1] = phi1;
        }
        else if (law_type == 3) // // Linear case + additional top value (discontinuity), 5 parameters expected: z0, z1, phi0, phi1, phi2
        // Line should be: 3 | 5 | z_0 | phi_0 | z_1 | phi_1 | phi_2
        // There are 7 records in the line
        {
            double z0 = cell_z_b[idx_cell] + records[idx_cell][2];
            double z1 = cell_z_b[idx_cell] + records[idx_cell][4];
            double phi0 = records[idx_cell][3];
            double phi1 = records[idx_cell][5];
            double phi2 = records[idx_cell][6];
            double dz = (z1 - z0) / (n_por_tab_max - 1);
            double dphi = (phi1 - phi0) / (n_por_tab_max - 1);
            cell_z_tab[idx_cell].resize(n_por_tab_max);
            phiW_tab[idx_cell].resize(n_por_tab_max);
            
            for (size_t z_level = 0; z_level < n_por_tab_max - 1; ++z_level)
            {
                cell_z_tab[idx_cell][z_level] = z0 + z_level * dz;
                phiW_tab[idx_cell][z_level] = phi0 + (z_level + 0.5) * dphi;
            }
            cell_z_tab[idx_cell][n_por_tab_max - 1] = z1;
            phiW_tab[idx_cell][n_por_tab_max - 1] = phi2;
        }
        //== Compute cell_z_tab & phiW_tab 

        //-- Compute thetaW_tab
        thetaW_tab[idx_cell] = {0.};
        double theta = 0;
        std::vector<double> z_tab = cell_z_tab[idx_cell];
        for (size_t z_level = 1; z_level < z_tab.size(); z_level++)
        {
            theta += (z_tab[z_level] - z_tab[z_level-1]) * phiW_tab[idx_cell][z_level-1];
            thetaW_tab[idx_cell].push_back(theta);
        }
        //== Compute thetaW_tab
        /*
        logInfo() << "Index, bottom, law type: " << i << cell_z_b[idx_cell] << records[idx_cell][0];
        logInfo() << "                       z:" << cell_z_tab[idx_cell];
        logInfo() << "                     phi:" << phiW_tab[idx_cell];
        logInfo() << "                   theta:" << thetaW_tab[idx_cell];
        */
    }
    return EXIT_SUCCESS;
}


bool swDdp::setPorosityEdgeBndMap(const std::vector<std::vector<double>>& records)
/**
 * @todo change into boo with EXIT_SUCCESS / FAILURE ?
 * @brief set porosity arrays (z, phi and theta) in boundary edges
 * @param[in] records : values of law_type, z and phi to be tabulated
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par modified attributes:
 *  - edge_bnd_z_tab
 *  - phiG_bnd_tab
 *  - thetaG_bnd_tab
 */
{
    for (size_t idx_edge_bnd = 0; idx_edge_bnd < edge_bnd_nb; idx_edge_bnd++)
    {
        //-- Compute edge_bnd_z_tab & phiG_bnd_tab
        double law_type = records[idx_edge_bnd][0];
        if (law_type == -2) // Arbitrary number n of couples (phi, z), these values are not re-tabulated
        // Line should be: law_type | n | z_0 | phi_0 | z_1 | phi_1 | ... | z_{n-1} | phi_{n-1}
        // There are 2n+2 records in the line
        {
            size_t n_couples = int(records[idx_edge_bnd][1]);
            logWarn() << "Using pre-tabulated porosities. Number of vertical porosity levels changed to" << n_couples;
            for (size_t z_level = 0; z_level < n_couples; z_level++)
            {
                double z = edge_bnd_midpoint_z[idx_edge_bnd] + records[idx_edge_bnd][2 + 2*z_level];
                double phi = records[idx_edge_bnd][3 + 2*z_level];
                edge_bnd_z_tab[idx_edge_bnd].push_back(z);
                phiG_bnd_tab[idx_edge_bnd].push_back(phi);
            }
        }
        else if (law_type == -1) // Arbitrary number n of couples (phi, z), these values are re-tabulated with n_por_tab_max levels
        // Line should be: law_type | n | z_0 | phi_0 | z_1 | phi_1 | ... | z_{n-1} | phi_{n-1}
        // There are 2n+2 records in the line
        {
            std::vector<double> phi_read = {};
            std::vector<double> z_read = {};
            size_t n_levels = int(records[idx_edge_bnd][1]);
            for (size_t z_level = 0; z_level < n_levels; z_level++)
            {
                z_read.push_back(edge_bnd_midpoint_z[idx_edge_bnd] + records[idx_edge_bnd][2 + 2*z_level]);
                phi_read.push_back(records[idx_edge_bnd][3 + 2*z_level]);
            }
            if (tabulatePorosity(n_por_tab_max, z_read, phi_read, edge_bnd_z_tab[idx_edge_bnd], phiG_bnd_tab[idx_edge_bnd]) != EXIT_SUCCESS)
            {
                logError() << "Error while tabulating porosity values with law type -1";
                return EXIT_FAILURE;
            }
        }
        else if (law_type == 0) // Arbitrary number n of (increasing) z levels, porosity to be interpolated between z_0 and z_{n-1}
        // Line should be: 0 | n | z_0 | z_1 | ... | z_{n-1}
        // There are n+2 records in the line
        {
            std::vector<double> phi_read = {};
            std::vector<double> z_read = {};
            size_t n_levels = int(records[idx_edge_bnd][1]);
            for (size_t z_level = 0; z_level < n_levels; z_level++)
            {
                z_read.push_back(edge_bnd_midpoint_z[idx_edge_bnd] + records[idx_edge_bnd][2 + z_level]);
                phi_read.push_back(z_level * 1.0 / (n_levels-1));
            }
            if (tabulatePorosity(n_por_tab_max, z_read, phi_read, edge_bnd_z_tab[idx_edge_bnd], phiG_bnd_tab[idx_edge_bnd]) != EXIT_SUCCESS)
            {
                logError() << "Error while tabulating porosity values with law type 0";
                return EXIT_FAILURE;
            }
        }
        else if (law_type == 1) // Uniform case, 2 parameters expected: z_0 and phi
        // Line should be: 1 | 2 | z_0 | phi
        // There are 4 records in the line
        {
            double z = edge_bnd_midpoint_z[idx_edge_bnd] + records[idx_edge_bnd][2];
            double phi = records[idx_edge_bnd][3];
            std::vector<double> phi_read = {};
            std::vector<double> z_read = {};
            z_read.push_back(z);
            z_read.push_back(z + 1.0);
            phi_read.push_back(phi);
            phi_read.push_back(phi);
            if (tabulatePorosity(n_por_tab_max, z_read, phi_read, edge_bnd_z_tab[idx_edge_bnd], phiG_bnd_tab[idx_edge_bnd]) != EXIT_SUCCESS)
            {
                logError() << "Error while tabulating porosity values with law type 1";
                return EXIT_FAILURE;
            }
            phiG_bnd_tab[idx_edge_bnd][0] = phi;
        }
        else if (law_type == 2) // Linear case, 4 parameters expected: z0, z1, phi0, phi1
        // Line should be: 2 | 4 | z_0 | phi_0 | z_1 | phi_1
        // There are 6 records in the line
        {
            double z0 = edge_bnd_midpoint_z[idx_edge_bnd] + records[idx_edge_bnd][2];
            double z1 = edge_bnd_midpoint_z[idx_edge_bnd] + records[idx_edge_bnd][4];
            double phi0 = records[idx_edge_bnd][3];
            double phi1 = records[idx_edge_bnd][5];
            double dz = (z1 - z0) / (n_por_tab_max - 1);
            double dphi = (phi1 - phi0) / (n_por_tab_max - 1);
            edge_bnd_z_tab[idx_edge_bnd].resize(n_por_tab_max);
            phiG_bnd_tab[idx_edge_bnd].resize(n_por_tab_max);
            
            for (size_t z_level = 0; z_level < n_por_tab_max - 1; ++z_level)
            {
                edge_bnd_z_tab[idx_edge_bnd][z_level] = z0 + z_level * dz;
                phiG_bnd_tab[idx_edge_bnd][z_level] = phi0 + (z_level + 0.5) * dphi;
            }
            edge_bnd_z_tab[idx_edge_bnd][n_por_tab_max - 1] = z1;
            phiG_bnd_tab[idx_edge_bnd][n_por_tab_max - 1] = phi1;
        }
        else if (law_type == 3) // // Linear case + additional top value (discontinuity), 5 parameters expected: z0, z1, phi0, phi1, phi2
        // Line should be: 3 | 5 | z_0 | phi_0 | z_1 | phi_1 | phi_2
        // There are 7 records in the line
        {
            double z0 = edge_bnd_midpoint_z[idx_edge_bnd] + records[idx_edge_bnd][2];
            double z1 = edge_bnd_midpoint_z[idx_edge_bnd] + records[idx_edge_bnd][4];
            double phi0 = records[idx_edge_bnd][3];
            double phi1 = records[idx_edge_bnd][5];
            double phi2 = records[idx_edge_bnd][6];
            double dz = (z1 - z0) / (n_por_tab_max - 1);
            double dphi = (phi1 - phi0) / (n_por_tab_max - 1);
            edge_bnd_z_tab[idx_edge_bnd].resize(n_por_tab_max);
            phiG_bnd_tab[idx_edge_bnd].resize(n_por_tab_max);
            
            for (size_t z_level = 0; z_level < n_por_tab_max - 1; ++z_level)
            {
                edge_bnd_z_tab[idx_edge_bnd][z_level] = z0 + z_level * dz;
                phiG_bnd_tab[idx_edge_bnd][z_level] = phi0 + (z_level + 0.5) * dphi;
            }
            edge_bnd_z_tab[idx_edge_bnd][n_por_tab_max - 1] = z1;
            phiG_bnd_tab[idx_edge_bnd][n_por_tab_max - 1] = phi2;
        }
        //== Compute edge_bnd_z_tab & phiG_bnd_tab 

        //-- Compute thetaG_bnd_tab
        thetaG_bnd_tab[idx_edge_bnd] = {0.};
        double theta = 0;
        std::vector<double> z_tab = edge_bnd_z_tab[idx_edge_bnd];
        for (size_t z_level = 1; z_level < z_tab.size(); z_level++)
        {
            theta += (z_tab[z_level] - z_tab[z_level-1]) * phiG_bnd_tab[idx_edge_bnd][z_level-1];
            thetaG_bnd_tab[idx_edge_bnd].push_back(theta);
        }
        //== Compute thetaG_bnd_tab
        /*
        logInfo() << "Edge_bnd index, bottom, law type: " << i << edge_bnd_midpoint_z[idx_edge_bnd] << records[idx_edge_bnd][0];
        logInfo() << "                               z:" << edge_bnd_z_tab[idx_edge_bnd];
        logInfo() << "                             phi:" << phiG_bnd_tab[idx_edge_bnd];
        logInfo() << "                           theta:" << thetaG_bnd_tab[idx_edge_bnd];
        */
    }
    return EXIT_SUCCESS;
}


bool swDdp::setPorosityEdgeIntMap(const std::vector<std::vector<double>>& records)
/**
 * @todo change into boo with EXIT_SUCCESS / FAILURE ?
 * @brief set porosity arrays (z, phi and theta) in internal edges
 * @param[in] records : values of law_type, z and phi to be tabulated
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par modified attributes:
 *  - edge_int_z_tab
 *  - phiG_int_tab
 *  - thetaG_int_tab
 *  - thetaG_int_hmin
 */
{
    for (size_t idx_edge_int = 0; idx_edge_int < edge_int_nb; idx_edge_int++)
    {
        //-- Compute edge_int_z_tab & phiG_int_tab
        double law_type = records[idx_edge_int][0];
        if (law_type == -2) // Arbitrary number n of couples (phi, z), these values are not re-tabulated
        // Line should be: law_type | n | z_0 | phi_0 | z_1 | phi_1 | ... | z_{n-1} | phi_{n-1}
        // There are 2n+2 records in the line
        {
            size_t n_couples = int(records[idx_edge_int][1]);
            logWarn() << "Using pre-tabulated porosities. Number of vertical porosity levels changed to" << n_couples;
            for (size_t z_level = 0; z_level < n_couples; z_level++)
            {
                double z = edge_int_midpoint_z[idx_edge_int] + records[idx_edge_int][2 + 2*z_level];
                double phi = records[idx_edge_int][3 + 2*z_level];
                edge_int_z_tab[idx_edge_int].push_back(z);
                phiG_int_tab[idx_edge_int].push_back(phi);
            }
        }
        else if (law_type == -1) // Arbitrary number n of couples (phi, z), these values are re-tabulated with n_por_tab_max levels
        // Line should be: law_type | n | z_0 | phi_0 | z_1 | phi_1 | ... | z_{n-1} | phi_{n-1}
        // There are 2n+2 records in the line
        {
            std::vector<double> phi_read = {};
            std::vector<double> z_read = {};
            phi_read = {};
            z_read = {}; 
            size_t n_levels = int(records[idx_edge_int][1]);
            for (size_t z_level = 0; z_level < n_levels; z_level++)
            {
                z_read.push_back(edge_int_midpoint_z[idx_edge_int] + records[idx_edge_int][2 + 2*z_level]);
                phi_read.push_back(records[idx_edge_int][3 + 2*z_level]);
            }
            if (tabulatePorosity(n_por_tab_max, z_read, phi_read, edge_int_z_tab[idx_edge_int], phiG_int_tab[idx_edge_int]) != EXIT_SUCCESS)
            {
                logError() << "Error while tabulating porosity values with law type -1";
                return EXIT_FAILURE;
            }
        }
        else if (law_type == 0) // Arbitrary number n of (increasing) z levels, porosity to be interpolated between z_0 and z_{n-1}
        // Line should be: 0 | n | z_0 | z_1 | ... | z_{n-1}
        // There are n+2 records in the line
        {
            std::vector<double> phi_read = {};
            std::vector<double> z_read = {};
            phi_read = {};
            z_read = {}; 
            size_t n_levels = int(records[idx_edge_int][1]);
            for (size_t z_level = 0; z_level < n_levels; z_level++)
            {
                z_read.push_back(edge_int_midpoint_z[idx_edge_int] + records[idx_edge_int][2 + z_level]);
                phi_read.push_back(z_level * 1.0 / (n_levels-1));
            }
            if (tabulatePorosity(n_por_tab_max, z_read, phi_read, edge_int_z_tab[idx_edge_int], phiG_int_tab[idx_edge_int]) != EXIT_SUCCESS)
            {
                logError() << "Error while tabulating porosity values with law type 0";
                return EXIT_FAILURE;
            }
        }
        else if (law_type == 1) // Uniform case, 2 parameters expected: z_0 and phi
        // Line should be: 1 | 2 | z_0 | phi
        // There are 4 records in the line
        {
            double z = edge_int_midpoint_z[idx_edge_int] + records[idx_edge_int][2];
            double phi = records[idx_edge_int][3];
            std::vector<double> phi_read = {};
            std::vector<double> z_read = {};
            z_read.push_back(z);
            z_read.push_back(z + 1.0);
            phi_read.push_back(phi);
            phi_read.push_back(phi);
            if (tabulatePorosity(n_por_tab_max, z_read, phi_read, edge_int_z_tab[idx_edge_int], phiG_int_tab[idx_edge_int]) != EXIT_SUCCESS)
            {
                logError() << "Error while tabulating porosity values with law type 1";
                return EXIT_FAILURE;
            }
            phiG_int_tab[idx_edge_int][0] = phi;
        }
        else if (law_type == 2) // Linear case, 4 parameters expected: z0, z1, phi0, phi1
        // Line should be: 2 | 4 | z_0 | phi_0 | z_1 | phi_1
        // There are 6 records in the line
        {
            double z0 = edge_int_midpoint_z[idx_edge_int] + records[idx_edge_int][2];
            double z1 = edge_int_midpoint_z[idx_edge_int] + records[idx_edge_int][4];
            double phi0 = records[idx_edge_int][3];
            double phi1 = records[idx_edge_int][5];
            double dz = (z1 - z0) / (n_por_tab_max - 1);
            double dphi = (phi1 - phi0) / (n_por_tab_max - 1);
            edge_int_z_tab[idx_edge_int].resize(n_por_tab_max);
            phiG_int_tab[idx_edge_int].resize(n_por_tab_max);
            
            for (size_t z_level = 0; z_level < n_por_tab_max - 1; ++z_level)
            {
                edge_int_z_tab[idx_edge_int][z_level] = z0 + z_level * dz;
                phiG_int_tab[idx_edge_int][z_level] = phi0 + (z_level + 0.5) * dphi;
            }
            edge_int_z_tab[idx_edge_int][n_por_tab_max - 1] = z1;
            phiG_int_tab[idx_edge_int][n_por_tab_max - 1] = phi1;
        }
        else if (law_type == 3) // // Linear case + additional top value (discontinuity), 5 parameters expected: z0, z1, phi0, phi1, phi2
        // Line should be: 3 | 5 | z_0 | phi_0 | z_1 | phi_1 | phi_2
        // There are 7 records in the line
        {
            double z0 = edge_int_midpoint_z[idx_edge_int] + records[idx_edge_int][2];
            double z1 = edge_int_midpoint_z[idx_edge_int] + records[idx_edge_int][4];
            double phi0 = records[idx_edge_int][3];
            double phi1 = records[idx_edge_int][5];
            double phi2 = records[idx_edge_int][6];
            double dz = (z1 - z0) / (n_por_tab_max - 1);
            double dphi = (phi1 - phi0) / (n_por_tab_max - 1);
            edge_int_z_tab[idx_edge_int].resize(n_por_tab_max);
            phiG_int_tab[idx_edge_int].resize(n_por_tab_max);
            for (size_t z_level = 0; z_level < n_por_tab_max - 1; ++z_level)
            {
                edge_int_z_tab[idx_edge_int][z_level] = z0 + z_level * dz;
                phiG_int_tab[idx_edge_int][z_level] = phi0 + (z_level + 0.5) * dphi;
            }
            edge_int_z_tab[idx_edge_int][n_por_tab_max - 1] = z1;
            phiG_int_tab[idx_edge_int][n_por_tab_max - 1] = phi2;
        }
        //== Compute edge_int_z_tab & phiG_int_tab 

        //-- Compute thetaG_int_tab
        thetaG_int_tab[idx_edge_int] = {0.};
        double theta = 0;
        std::vector<double> z_tab = edge_int_z_tab[idx_edge_int];
        for (size_t z_level = 1; z_level < z_tab.size(); z_level++)
        {
            theta += (z_tab[z_level] - z_tab[z_level-1]) * phiG_int_tab[idx_edge_int][z_level-1];
            thetaG_int_tab[idx_edge_int].push_back(theta);
        }
        //== Compute thetaG_int_tab
        //-- Compute thetaG_int_hmin
        thetaG_int_hmin[idx_edge_int] = thetaGIntOfZ(idx_edge_int, z_tab[0] + h_min_ddp);
        //== Compute thetaG_int_hmin
    }
    return EXIT_SUCCESS;
}


bool swDdp::setPorosityParameters(void)
/**
 * @brief Reads the distribution of DDP parameters from file and sets the corresponding tables
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par modified attributes:
 *  - cell_z_tab, edge_bnd_tab, edge_int_tab
 *  - phiW_tab, phiG_bnd_tab, phiG_bnd_tab
 *  - thetaW_tab, thetaG_bnd_tab, thetaG_int_tab
 */
{
    logInfo() << "swDdp::setPorosityParameters - 1/3 - reading cell porosity map";
    ///- Cell porosity: Set up pre-defined input file name as defined in file.h
    std::string file_name;
    file_name = swFile::input_ddp_porosity_cell_file;
    file_name = m_io->absolutePath(file_name);
    std::vector<std::string> headers;
    std::vector<std::vector<double>> records;
    if (readPorosityFile(file_name, cell_nb, headers, records) != EXIT_SUCCESS)
    {
        logError() << "Error while reading porosity file";
        return EXIT_FAILURE;
    }
    cell_z_tab.resize(cell_nb);
    phiW_tab.resize(cell_nb);
    thetaW_tab.resize(cell_nb);
    if (setPorosityCellMap(records) != EXIT_SUCCESS)
    {
        logError() << "Error while setting cell porosity map";
        return EXIT_FAILURE;
    }
    logInfo() << "                             - 1/3 done";
    //== Cell porosity: Set up pre-defined input file name as defined in file.h

    ///- Edge_bnd porosity: Set up pre-defined input file name as defined in file.h
    logInfo() << "swDdp::setPorosityParameters - 2/3 - reading edge_bnd porosity map";
    file_name = swFile::input_ddp_porosity_edge_bnd_file;
    file_name = m_io->absolutePath(file_name);
    if (readPorosityFile(file_name, edge_bnd_nb, headers, records) != EXIT_SUCCESS)
    {
        logError() << "Error while reading porosity file";
        return EXIT_FAILURE;
    }
    edge_bnd_z_tab.resize(edge_bnd_nb);
    phiG_bnd_tab.resize(edge_bnd_nb);
    thetaG_bnd_tab.resize(edge_bnd_nb);
    if (setPorosityEdgeBndMap(records) != EXIT_SUCCESS)
    {
        logError() << "Error while setting boundary edges porosity map";
        return EXIT_FAILURE;
    }
    logInfo() << "                             - 2/3 done";
    //== Edge_bnd porosity: Set up pre-defined input file name as defined in file.h

    ///- Edge_int porosity: Set up pre-defined input file name as defined in file.h
    logInfo() << "swDdp::setPorosityParameters - 3/3 - reading edge_int porosity map";
    file_name = swFile::input_ddp_porosity_edge_int_file;
    file_name = m_io->absolutePath(file_name);
    if (readPorosityFile(file_name, edge_int_nb, headers, records) != EXIT_SUCCESS)
    {
        logError() << "Error while reading porosity file";
        return EXIT_FAILURE;
    }
    edge_int_z_tab.resize(edge_int_nb);
    phiG_int_tab.resize(edge_int_nb);
    thetaG_int_tab.resize(edge_int_nb);
    thetaG_int_hmin.resize(edge_int_nb);
    if (setPorosityEdgeIntMap(records) != EXIT_SUCCESS)
    {
        logError() << "Error while setting internal edges porosity map";
        return EXIT_FAILURE;
    }
    logInfo() << "                             - 3/3 done";
    //== Edge_int porosity: Set up pre-defined input file name as defined in file.h
    return EXIT_SUCCESS;
}


/** @todo: adapt code */
void swDdp::setRunoffParameters(void)
/**
 * @brief read information from runoff file (could be full rainfall, infiltration coeff, runoff coeff, SCS model)
 * @par modified attributes:
 *  - runoff_param
 *  - scs_i_a
 *  - scs_s
 */
{
    logInfo() << "\n \n";
    logInfo() << "swDdp::init - reading runoff parameter map";
    std::string file_name;
    std::vector<std::string> variable_name;
    std::vector<std::vector<double>> param_array;
    if (model_runoff == "full_rainfall")
    {
        // No file import required
    }
    else if (model_runoff == "infiltration_coefficient")
    {
        file_name = m_io->absolutePath(swFile::input_runoff_infiltration_coefficient_map);
        runoff_param.resize(cell_nb);
        swIO::readVectorParameterFile(file_name, 1, cell_nb, variable_name, param_array);
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            runoff_param[idx_cell] = param_array[0][idx_cell];
        }
    }
    else if (model_runoff == "runoff_coefficient")
    {
        file_name = m_io->absolutePath(swFile::input_runoff_runoff_coefficient_map);
        runoff_param.resize(cell_nb);
        swIO::readVectorParameterFile(file_name, 1, cell_nb, variable_name, param_array);
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            runoff_param[idx_cell] = param_array[0][idx_cell];
        }
    }
    else if (model_runoff == "scs") 
    {
        //-- Read SCS model parameters
        file_name = m_io->absolutePath(swFile::input_runoff_scs_parameter_map);
        swIO::readVectorParameterFile(file_name, 2, cell_nb, variable_name, param_array);
        scs_i_a.resize(cell_nb);
        scs_s.resize(cell_nb);
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            scs_i_a[idx_cell] = param_array[0][idx_cell];
            scs_s[idx_cell] = param_array[1][idx_cell];
        }
        //== Read SCS model parameters
        //-- Resize variables
        scs_h_a.resize(cell_nb);
        scs_h_s.resize(cell_nb);
        //== Resize variables
    }
    else
    {
        logError() << "Runoff model not handled = " << model_runoff;
        logFatal(__FILE__, __func__, __LINE__);
    }
    logInfo() << "             - Runoff parameter map read \n";
}


/** @todo: adapt code */
void swDdp::setSingularHeadLossParameters(void)
/**
 * @brief read information from singular head loss file
 * @par modified attributes:
 *  - singular_head_loss_h_sing
 *  - singular_head_loss_coeff_direction_1
 *  - singular_head_loss_coeff_direction_2
 *  - singular_head_loss_coeff_alpha
 */
{
    logInfo() << "swDdp::init - reading singular head loss map";
    std::string file_name = m_io->absolutePath(swFile::input_singular_head_loss_map);
    std::vector<std::string> variable_name;
    std::vector<std::vector<double>> param_array;
    swIO::readVectorParameterFile(file_name, 4, cell_nb, variable_name, param_array);
    singular_head_loss_h_sing.resize(cell_nb);
    singular_head_loss_coeff_direction_1.resize(cell_nb);
    singular_head_loss_coeff_direction_2.resize(cell_nb);
    singular_head_loss_alpha.resize(cell_nb);
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        singular_head_loss_h_sing[idx_cell] = param_array[0][idx_cell];
        singular_head_loss_coeff_direction_1[idx_cell] = param_array[1][idx_cell];
        singular_head_loss_coeff_direction_2[idx_cell] = param_array[2][idx_cell];
        singular_head_loss_alpha[idx_cell] = param_array[3][idx_cell] * atan(1.0) / 45.0;  // Convert degrees to radians
    }
    logInfo() << "             - friction map read \n";
}


/** @todo: adapt code */
void swDdp::setTransportDiffusionParameters(void)
/**
 * @brief sets transport diffusion array (cellwise)
 * @par modified attributes:
 *  - transport_scalar_diffusion_parameter
 */
{
}


double swDdp::specificForce(const double &theta_velocity, const double &velocity, const double &pressure_force)
/**
 * @brief Compute specific force
 * @param theta_velocity 
 * @param velocity 
 * @param pressure 
 * @return  
 */
{
    return theta_velocity * velocity + pressure_force;
}


/** @todo: adapt code */
void swDdp::storeResultMaps(const std::string& map_storage_mode, double const& time)
/**
 * @brief Writes the result map file
 * @param[in] map_storage_mode: string to define the output verbosity (default or all)
 * @param[in] time : current time
 * @par Modified attributes:
 * - content of output files
 */
{
    if (map_storage_mode != "none")
    {
        ///-- Set up file name
        std::string file_name =swFile::output_hydro_map_file_root;
        m_io->computeFileNameFromTime(t, file_name);
        //== Set up file name

        ///-- Fill in vectors to be written
        std::vector<std::string> param_name = {};
        if (map_storage_mode == "default")
        {
            param_name = {"x","y","zb", 
                          "water_depth", "unit_discharge_x", "unit_discharge_y ", "unit_disch_norm  "};
        }
        else if (map_storage_mode == "all")
        {
            param_name = {"x", "y", "zb", 
                          "water_depth", "unit_discharge_x", "unit_discharge_y", "unit_disch_norm",
                          "surface_elevation", "velocity_x", "velocity_y", "velocity_norm",
                          "froude", "specific_head", "total_head",
                          "theta", "theta_velocity_x", "theta_velocity_y", "theta_veloc_norm"};
        }
        std::vector<std::vector<double>> param_array;
        param_array.resize(param_name.size());


        //-- Parts of the table that will always be filled in 
        for (size_t n = 0; n < cell_nb; n++) // Instructions are faster than copy
        {
            {
                param_array[0].push_back(cell_centre_x[n]);
                param_array[1].push_back(cell_centre_y[n]);
                param_array[2].push_back(cell_z_tab[n][0]);
                double theta = thetaW[n];
                double z = zOfThetaW(n, theta);
                double h = std::max(0.0, z - cell_z_tab[n][0]);
                param_array[3].push_back(h);
                double u = 0.0;
                double v = 0.0;
                if (theta > small_real)
                {
                    u = thetaW_velocity_x[n] / theta;
                    v = thetaW_velocity_y[n] / theta;
                }
                param_array[4].push_back(u * h);
                param_array[5].push_back(v * h);
                double q_norm = h * std::sqrt(u * u + v * v);
                param_array[6].push_back(q_norm);
            }
            //== Parts of the table that will be filled in "default", "full" and "nonzero" modes
        }

        //-- Optional : only in 'all' mode
        if (map_storage_mode == "all"  )
        {
            for (size_t n = 0; n < cell_nb; n++) // Instructions are faster than copy
            {
                double theta = thetaW[n];
                double z = zOfThetaW(n, theta);
                double h = std::max(0.0, z - cell_z_tab[n][0]);
                param_array[7].push_back(z);
                double u = 0.0;
                double v = 0.0;
                if (theta > small_real)
                {
                    u = thetaW_velocity_x[n] / theta;
                    v = thetaW_velocity_y[n] / theta;
                }
                param_array[8].push_back(u);
                param_array[9].push_back(v);
                double v_norm =  std::sqrt(u * u + v * v);
                param_array[10].push_back(v_norm);
                double phi = phiWOfZ(n, z);
                double c = pressureWaveCelerity(theta, phi);
                param_array[11].push_back(v_norm / c);
                param_array[12].push_back(h + 0.5 * v_norm * v_norm / gravity);
                param_array[13].push_back(z + 0.5 * v_norm * v_norm / gravity);
                param_array[14].push_back(theta);
                param_array[15].push_back(u * theta);
                param_array[16].push_back(v * theta);
                param_array[17].push_back(v_norm * theta);
            }
        }
        //== Optional : only in 'all' mode

        //== Fill in vectors to be written

        ///-- Call writer
        if ((output_format == "txt") || (output_format == "all"))
        {
            if (swIO::createTxtFileMapDecimal(file_name, param_name, param_array, n_decimal) != EXIT_SUCCESS)
            {
                logError() << "swGeneric::writeResultMap - Exit Failure in swIO::createTxtMapDecimal \n";
                logFatal(__FILE__, __func__, __LINE__);
            }
        }
        if ((output_format == "bin") || (output_format == "all"))
        {
            std::string filename = m_io->absolutePath(swFile::output_map_binary_file);
            if (swIO::fillBinFileMapDecimal(filename, t, param_array) != EXIT_SUCCESS)
            {
                logError() << "swGeneric::writeResultMap - Exit Failure in swIO::fillBinFileMapDecimal \n";
                logFatal(__FILE__, __func__, __LINE__);
            }
        }
        //== Call writer
    }
}


/** @todo: adapt code */
void swDdp::storeResultProbes(double const& time)
/**
 * @brief Stores the probes values on one single line
 * @param[in] time : current time
 * @par Modified attributes:
 * - content of output files
 */
{   
    /// Line starts with time
    std::string probes_file_line = swString::numberToStr(time, 'e', 10) + "\t";
    
    ///-- For each probe, print z_b, h, q and r
    for (size_t p = 0 ; p < probes_nb ; p++)
    {
        size_t id_cell = probe_cell_idx[p];
        double theta = thetaW[id_cell];
        double z = zOfThetaW(id_cell, theta);
        double h = std::max(0.0, z - cell_z_tab[id_cell][0]);
        double u = thetaW_velocity_x[id_cell] / std::max(small_real, theta);
        double v = thetaW_velocity_y[id_cell] / std::max(small_real, theta);
        probes_file_line += swString::numberToStr(cell_z_tab[id_cell][0], 'e', 10) + "\t" + swString::numberToStr(h, 'e', 10) 
            + "\t" + swString::numberToStr(h * u, 'e', 10) + "\t" + swString::numberToStr(h * v, 'e', 10) + "\t";
    }
    probes_file_line += "\n";
    //== For each probe, print z_b, h, q and r
        
    ///-- Call writer to print line in probe file
    if (swIO::writeTxtLine(probes_output_stream, probes_file_line) != EXIT_SUCCESS)
    {
        logError() << "swDdp::storeResultProbes - Exit Failure in swIO::writeTxtLine \n";
        logFatal(__FILE__, __func__, __LINE__);
    }
    //== Call writer
}


bool swDdp::tabulatePorosity(size_t ntab, const std::vector<double> &z, const std::vector<double> &phi, 
                                     std::vector<double>& z_tab, std::vector<double>& phi_tab)
/**
 * @brief Tabulates z and phi to final arrays 
 * @param[in] ntab : number of levels in the output arrays
 * @param[in] z : input array for altitudes (same size as phi)
 * @param[in] phi : input array for porosities (same size as z)
 * @param[out] z_tab : output tabulated values for altitudes (size = ntab)
 * @param[out] phi_tab : output tabulated values for altitudes (size = ntab)
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 */
{
    size_t nb_values = z.size();
    phi_tab.resize(ntab);
    z_tab.resize(ntab);
    //-- interpolate phi
    double dz = (z[nb_values - 1] - z[0]) / (ntab - 1);
    for (size_t level = 0; level < ntab; ++level)
    {
        z_tab[level] = z[0] + level * dz;
    }

    std::vector<double> phi_interp;
    phi_interp.resize(ntab);
    phi_interp[0] = 0;
    
    size_t index = 0;
    for (size_t level = 1; level < ntab; ++level)
    {
        while (z_tab[level] > z[index])
        {
            ++index;
        }
        phi_interp[level] = phi[index - 1] +
                (phi[index] - phi[index - 1]) *
                (z_tab[level] - z[index - 1]) /
                (z[index] - z[index - 1]);
    }
    //== interpolate phi
    
    //-- compute phi
    for (size_t level = 0; level < ntab - 1; ++level)
    {
        phi_tab[level] = (phi_interp[level + 1] + phi_interp[level]) * 0.5;
    }
    phi_tab[ntab - 1] = phi[nb_values - 1];
    //== compute phi

    return EXIT_SUCCESS;
}


double swDdp::thetaGBndOfZ(const size_t& idx_edge_bnd, const double& z)
/**
 * @brief Compute theta(z) in any given boundary edge, where theta is pre-tabulated  from the porosity data
 * @param idx_edge_bnd: boundary edge index in which the integral porosity is computed
 * @param z: altitude where the boundary edge integral porosity is evaluated
 * @return evaluated integral porosity in boundary edge idx_edge_bnd, at altitude z
 */
{
    std::vector<double> z_tab = edge_bnd_z_tab[idx_edge_bnd];
    if (z <= z_tab[0])
    {
        return 0.;
    }
    else
    {
        size_t level = 0;
        while ((level < z_tab.size()) && (z > z_tab[level]))
        {
            ++level;
        }
        double theta = thetaG_bnd_tab[idx_edge_bnd][level - 1];
        if (z > z_tab[level - 1])
        {
            theta += (z - z_tab[level - 1]) * phiG_bnd_tab[idx_edge_bnd][level - 1];
        }
        return theta;
    }
}


double swDdp::thetaGIntOfZ(const size_t& idx_edge_int, const double& z)
/**
 * @brief Compute theta(z) in any given internal edge, where theta is pre-tabulated from the porosity data
 * @param idx_edge_int: internal edge index in which the integral porosity is computed
 * @param z: altitude where the internal edge integral porosity is evaluated
 * @return evaluated integral porosity in internal edge idx_edge_int, at altitude z
 */
{
    std::vector<double> z_tab = edge_int_z_tab[idx_edge_int];
    if (z <= z_tab[0])
    {
        return 0.;
    }
    else
    {
        size_t level = 0;
        while ((level < z_tab.size()) && (z > z_tab[level]))
        {
            ++level;
        }
        double theta = thetaG_int_tab[idx_edge_int][level - 1];
        if (z > z_tab[level - 1])
        {
            theta += (z - z_tab[level - 1]) * phiG_int_tab[idx_edge_int][level - 1];
        }
        return theta;
    }
}


double swDdp::thetaWOfZ(const size_t& idx_cell, const double& z)
/**
 * @brief Compute theta(z) in any given cell, where theta is pre-tabulated from the porosity data
 * @param idx_cell: cell index in which the integral porosity is computed
 * @param z: altitude where the cell integral porosity is evaluated
 * @return evaluated integral porosity in cell idx_cell, at altitude z
 */
{
    std::vector<double> z_tab = cell_z_tab[idx_cell];
    if (z <= z_tab[0])
    {
        return 0.;
    }
    else
    {
        size_t level = 0;
        while ((level < z_tab.size()) && (z > z_tab[level]))
        {
            ++level;
        }
        double theta = thetaW_tab[idx_cell][level - 1];
        if (z > z_tab[level - 1])
        {
            theta += (z - z_tab[level - 1]) * phiW_tab[idx_cell][level - 1];
        }
        return theta;
    }
}


/** @todo: adapt code */
void swDdp::writeLogBoussinesq(void)
/**
 * @brief Writes down the boussinesq coeff log file
 * @return nothing
 * @par Modified attributes:
 * - parameter_names
 * - parameter_values
 * - content of log files
 */
{
    logInfo() << "swSp - writeLogBoussinesq";
    //-- Create the string vector that will be used for the header
    std::vector<std::string> parameter_names = {"x", "y", "idx_cell"};
    parameter_names.push_back("Boussinesq_coefficient");
    //== Create the string vector that will be used for the header


    //-- Create the string table for parameter values
    std::vector<std::vector<std::string>> parameter_values = {};
    parameter_values.resize(parameter_names.size());
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        parameter_values[0].push_back(std::to_string(cell_centre_x[idx_cell]));
        parameter_values[1].push_back(std::to_string(cell_centre_y[idx_cell]));
        parameter_values[2].push_back(std::to_string(idx_cell));
        parameter_values[3].push_back(std::to_string(boussinesq[idx_cell]));
    }
    //== Create the string table for parameter values

    //-- Create file and write values
    std::string file_name = "log/Boussinesq_log_file.txt";
    file_name = m_io->absolutePath(file_name);
    swIO::createFileMapString(file_name, parameter_names, parameter_values);
    //== Create file and write values
    logInfo() << "swSp - writeLogBoussinesq done";
}


/** @todo: adapt code */
void swDdp::writeLogEdgeFactors(void)
/**
 * @brief Writes log file for edge factors
 * @par Modified attributes:
 * - parameter_names
 * - parameter_values
 * - content of log files
 */
{
    logInfo() << "swSp - writeLogEdgeFactors";

    //-- Boundary edges
    //---- Create the string vector that will be used for the header
    std::vector<std::string> parameter_names = {"x", "y", "z","idx_edge"};
    parameter_names.push_back("fact_1_F1");
    parameter_names.push_back("fact_2_F2");
    parameter_names.push_back("fact_2_F3");
    parameter_names.push_back("fact_2_S2");
    parameter_names.push_back("fact_1_F1");
    parameter_names.push_back("fact_3_F3");
    parameter_names.push_back("fact_3_S2");
    //==== Create the string vector that will be used for the header


    //---- Create the string table for parameter values
    std::vector<std::vector<std::string>> parameter_values = {};
    parameter_values.resize(parameter_names.size());
    for (size_t idx_edge_bnd = 0; idx_edge_bnd < edge_bnd_nb; idx_edge_bnd++)
    {
        parameter_values[0].push_back(std::to_string(edge_bnd_midpoint_x[idx_edge_bnd]));
        parameter_values[1].push_back(std::to_string(edge_bnd_midpoint_y[idx_edge_bnd]));
        parameter_values[2].push_back(std::to_string(edge_bnd_midpoint_z[idx_edge_bnd]));
        parameter_values[3].push_back(std::to_string(idx_edge_bnd));
        parameter_values[4].push_back(std::to_string(edge_bnd_factor_1_F1[idx_edge_bnd]));
        parameter_values[5].push_back(std::to_string(edge_bnd_factor_2_F2[idx_edge_bnd]));
        parameter_values[6].push_back(std::to_string(edge_bnd_factor_2_F3[idx_edge_bnd]));
        parameter_values[7].push_back(std::to_string(edge_bnd_factor_2_S2[idx_edge_bnd]));
        parameter_values[8].push_back(std::to_string(edge_bnd_factor_3_F2[idx_edge_bnd]));
        parameter_values[9].push_back(std::to_string(edge_bnd_factor_3_F3[idx_edge_bnd]));
        parameter_values[10].push_back(std::to_string(edge_bnd_factor_3_S2[idx_edge_bnd]));
    }
    //==== Create the string table for parameter values

    //---- Create file and write values
    std::string file_name = swFile::log_edge_bnd_factors_file;
    file_name = m_io->absolutePath(file_name);
    swIO::createFileMapString(file_name, parameter_names, parameter_values);
    //==== Create file and write values
    //== Boundary edges

    //-- Internal edges
    //---- Create the string vector that will be used for the header
    parameter_names = {"x", "y", "z", "idx_edge_int"};
    parameter_names.push_back("fact_1_F1_L");
    parameter_names.push_back("fact_1_F1_R");
    parameter_names.push_back("fact_2_F2_L");
    parameter_names.push_back("fact_2_F2_R");
    parameter_names.push_back("fact_2_F3_L");
    parameter_names.push_back("fact_2_F3_R");
    parameter_names.push_back("fact_2_S2_L");
    parameter_names.push_back("fact_2_S2_R");
    parameter_names.push_back("fact_1_F1_L");
    parameter_names.push_back("fact_1_F1_R");
    parameter_names.push_back("fact_3_F3_L");
    parameter_names.push_back("fact_3_F3_R");
    parameter_names.push_back("fact_3_S2_L");
    parameter_names.push_back("fact_3_S2_R");
    //==== Create the string vector that will be used for the header


    //---- Create the string table for parameter values
    parameter_values = {};
    parameter_values.resize(parameter_names.size());
    for (size_t idx_edge_int = 0; idx_edge_int < edge_int_nb; idx_edge_int++)
    {
        parameter_values[0].push_back(std::to_string(edge_int_midpoint_x[idx_edge_int]));
        parameter_values[1].push_back(std::to_string(edge_int_midpoint_y[idx_edge_int]));
        parameter_values[2].push_back(std::to_string(edge_int_midpoint_z[idx_edge_int]));
        parameter_values[3].push_back(std::to_string(idx_edge_int));
        parameter_values[4].push_back(std::to_string(edge_int_factor_1_F1_L[idx_edge_int]));
        parameter_values[5].push_back(std::to_string(edge_int_factor_1_F1_R[idx_edge_int]));
        parameter_values[6].push_back(std::to_string(edge_int_factor_2_F2_L[idx_edge_int]));
        parameter_values[7].push_back(std::to_string(edge_int_factor_2_F2_R[idx_edge_int]));
        parameter_values[8].push_back(std::to_string(edge_int_factor_2_F3_L[idx_edge_int]));
        parameter_values[9].push_back(std::to_string(edge_int_factor_2_F3_R[idx_edge_int]));
        parameter_values[10].push_back(std::to_string(edge_int_factor_2_S2_L[idx_edge_int]));
        parameter_values[11].push_back(std::to_string(edge_int_factor_2_S2_R[idx_edge_int]));
        parameter_values[12].push_back(std::to_string(edge_int_factor_3_F2_L[idx_edge_int]));
        parameter_values[13].push_back(std::to_string(edge_int_factor_3_F2_R[idx_edge_int]));
        parameter_values[14].push_back(std::to_string(edge_int_factor_3_F3_L[idx_edge_int]));
        parameter_values[15].push_back(std::to_string(edge_int_factor_3_F3_R[idx_edge_int]));
        parameter_values[16].push_back(std::to_string(edge_int_factor_3_S2_L[idx_edge_int]));
        parameter_values[17].push_back(std::to_string(edge_int_factor_3_S2_R[idx_edge_int]));
    }
    //==== Create the string table for parameter values

    //---- Create file and write values
    file_name = swFile::log_edge_int_factors_file;
    file_name = m_io->absolutePath(file_name);
    swIO::createFileMapString(file_name, parameter_names, parameter_values);
    //==== Create file and write values
    //== Boundary edges

    logInfo() << "swSp - writeLogEdgeFactors done";
}


/** @todo: adapt code */
void swDdp::writeLogFriction(void)
/**
 * @brief Writes log file for friction
 * @par Modified attributes:
 *  - parameter_names
 *  - parameter_values
 *  - content of log files
 */
{
    logInfo() << "swSp - writeLogFriction - model_friction = " << model_friction;
    //-- Create the string vector that will be used for the header
    std::vector<std::string> parameter_names = {"x", "y", "idx_cell"};
    parameter_names.push_back(model_friction + "_coeff_1");
    parameter_names.push_back(model_friction + "_coeff_2");
    parameter_names.push_back(model_friction + "_alpha");
    //== Create the string vector that will be used for the header


    //-- Create the string table for parameter values
    std::vector<std::vector<std::string>> parameter_values = {};
    parameter_values.resize(parameter_names.size());
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        parameter_values[0].push_back(std::to_string(cell_centre_x[idx_cell]));
        parameter_values[1].push_back(std::to_string(cell_centre_y[idx_cell]));
        parameter_values[2].push_back(std::to_string(idx_cell));
        parameter_values[3].push_back(std::to_string(friction_coeff_direction_1[idx_cell]));
        parameter_values[4].push_back(std::to_string(friction_coeff_direction_2[idx_cell]));
        parameter_values[5].push_back(std::to_string(friction_alpha[idx_cell]));
    }
    //== Create the string table for parameter values

    //-- Create file and write values
    std::string file_name = swFile::log_friction_file;
    file_name = m_io->absolutePath(file_name);
    swIO::createFileMapString(file_name, parameter_names, parameter_values);
    //== Create file and write values
    logInfo() << "swSp - writeLogFriction done";
}


/** @todo: adapt code */
void swDdp::writeLogInfiltration(void)
/**
 * @brief Writes log file for infiltration
 * @par Modified attributes:
 *  - parameter_names
 *  - parameter_values
 *  - content of log files
 */
{
    logInfo() << "swDdp - writeLogInfiltration - model_infiltration = " << model_infiltration;
    //-- Create the string vector that will be used for the header
    std::vector<std::string> parameter_names = { "x", "y", "idx_cell" };
    parameter_names.push_back(model_infiltration + "_parameter");
    //== Create the string vector that will be used for the header


    //-- Create the string table for parameter values
    std::vector<std::vector<std::string>> parameter_values = {};
    parameter_values.resize(parameter_names.size());
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        parameter_values[0].push_back(std::to_string(cell_centre_x[idx_cell]));
        parameter_values[1].push_back(std::to_string(cell_centre_y[idx_cell]));
        parameter_values[2].push_back(std::to_string(idx_cell));
        parameter_values[3].push_back(std::to_string(infiltration_param[idx_cell]));
    }
    //== Create the string table for parameter values

    //-- Create file and write values
    std::string file_name = swFile::log_infiltration_file;
    file_name = m_io->absolutePath(file_name);
    swIO::createFileMapString(file_name, parameter_names, parameter_values);
    //== Create file and write values
    logInfo() << "swDdp - writeLogInfiltration done";
}


/** @todo: adapt code */
void swDdp::writeLogMomentumDiffusion(void)
/**
 * @brief Writes log file for runoff
 * @param
 * @todo Write the code
 */
{
    logInfo() << "Not implemented yet";
}


bool swDdp::writeLogPorosity(void)
/**
 * @brief Writes log file for DDP porosities
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 * @par Modified attributes:
 * - parameter_names
 * - parameter_values
 * - content of log files
 */
{
    // Cells
    std::ofstream file_w;
    std::string file_name = swFile::log_ddp_porosity_cell_file;
    file_name = m_io->absolutePath(file_name);
    file_w.open(file_name);
    if(!file_w)

    {
        logError() << "    Unable to open the DDP PhiW log file:" << swFile::log_ddp_porosity_cell_file;
        return EXIT_FAILURE;
    }
    logInfo() << "Writing DDP PhiW log file:" << swFile::log_ddp_porosity_cell_file;

    file_w << "x\ty\tNbTab\n";
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        double x = cell_centre_x[idx_cell];
        double y = cell_centre_y[idx_cell];
        size_t nb_tab = cell_z_tab[idx_cell].size();
        file_w << std::to_string(x) << "\t"
                    << std::to_string(y) << "\t"
                    << nb_tab << "\n";
        file_w << "\tzTab";
        for (size_t j = 0; j < nb_tab; ++j)
        {
            file_w << "\t" << std::to_string(cell_z_tab[idx_cell][j]);
        }
        file_w << "\n";
        file_w << "\tPhi ";
        for (size_t j = 0; j < nb_tab; ++j)
        {
            file_w << "\t" << std::to_string(phiW_tab[idx_cell][j]);
        }
        file_w << "\n";
        file_w << "\tTheta";
        for (size_t j = 0; j < nb_tab; ++j)
        {
            file_w << "\t" << std::to_string(thetaW_tab[idx_cell][j]);
        }
        file_w << "\n";
    }
    file_w.close();
    logEmpt() << "...done";

    // Edges_bnd
    file_name = swFile::log_ddp_porosity_edge_bnd_file;
    file_name = m_io->absolutePath(file_name);
    file_w.open(file_name);
    if(!file_w)
    {
        logError() << "    Unable to open the DDP PhiG BND log file:" << swFile::log_ddp_porosity_edge_bnd_file;
        return EXIT_FAILURE;
    }
    logInfo() << "Writing DDP PhiG BND log file:" << swFile::log_ddp_porosity_edge_bnd_file;

    file_w << "x\ty\tNbTab\tzTab\tPhi\tTheta\n";
    for (size_t idx_edge_bnd = 0; idx_edge_bnd < edge_bnd_nb; idx_edge_bnd++)
    {
        double x = edge_bnd_midpoint_x[idx_edge_bnd];
        double y = edge_bnd_midpoint_y[idx_edge_bnd];
        size_t nb_tab = edge_bnd_z_tab[idx_edge_bnd].size();
        file_w << std::to_string(x) << "\t"
                    << std::to_string(y) << "\t"
                    << nb_tab << "\n";
        file_w << "\tzTab";
        for (size_t j = 0; j < nb_tab; ++j)
        {
            file_w << "\t" << std::to_string(edge_bnd_z_tab[idx_edge_bnd][j]);
        }
        file_w << "\n";
        file_w << "\tPhi ";
        for (size_t j = 0; j < nb_tab; ++j)
        {
            file_w << "\t" << std::to_string(phiG_bnd_tab[idx_edge_bnd][j]);
        }
        file_w << "\n";
        file_w << "\tTheta";
        for (size_t j = 0; j < nb_tab; ++j)
        {
            file_w << "\t" << std::to_string(thetaG_bnd_tab[idx_edge_bnd][j]);
        }
        file_w << "\n";
    }
    file_w.close();
    logEmpt() << "...done";

    // Edges_int
    file_name = swFile::log_ddp_porosity_edge_int_file;
    file_name = m_io->absolutePath(file_name);
    file_w.open(file_name);
    if(!file_w)
    {
        logError() << "    Unable to open the DDP PhiG INT log file:" << swFile::log_ddp_porosity_edge_int_file;
        return EXIT_FAILURE;
    }
    logInfo() << "Writing DDP PhiG INT log file:" << swFile::log_ddp_porosity_edge_int_file;

    file_w << "x\ty\tNbTab\tzTab\tPhi\tTheta\n";
    for (size_t idx_edge_int = 0; idx_edge_int < edge_int_nb; idx_edge_int++)
    {
        double x = edge_int_midpoint_x[idx_edge_int];
        double y = edge_int_midpoint_y[idx_edge_int];
        size_t nb_tab = edge_int_z_tab[idx_edge_int].size();
        file_w << std::to_string(x) << "\t"
                    << std::to_string(y) << "\t"
                    << nb_tab << "\n";
        file_w << "\tzTab";
        for (size_t j = 0; j < nb_tab; ++j)
        {
            file_w << "\t" << std::to_string(edge_int_z_tab[idx_edge_int][j]);
        }
        file_w << "\n";
        file_w << "\tPhi ";
        for (size_t j = 0; j < nb_tab; ++j)
        {
            file_w << "\t" << std::to_string(phiG_int_tab[idx_edge_int][j]);
        }
        file_w << "\n";
        file_w << "\tTheta";
        for (size_t j = 0; j < nb_tab; ++j)
        {
            file_w << "\t" << std::to_string(thetaG_int_tab[idx_edge_int][j]);
        }
        file_w << "\n";
    }
    file_w.close();
    logEmpt() << "...done";
    return EXIT_SUCCESS;
}


/** @todo: adapt code */
void swDdp::writeLogRunoff(void)
/**
 * @brief Writes log file for runoff
 * @par Modified attributes:
 *  - parameter_names
 *  - parameter_values
 *  - content of log files
 */
{
    logInfo() << "swSp - writeLogRunoff - model_runoff = " << model_runoff;
    //-- Create the string vector that will be used for the header
    std::vector<std::string> parameter_names = { "x", "y", "idx_cell", "precip_station"};
    if ((model_runoff == "infiltration_coefficient") ||
        (model_runoff == "runoff_coefficient" )         )
    {
        parameter_names.push_back(model_runoff + "_parameter");
    }
    else if (model_runoff == "scs")
    {
        parameter_names.push_back("scs_I_a");
        parameter_names.push_back("scs_S");
    }
    //== Create the string vector that will be used for the header


    //-- Create the string table for parameter values
    std::vector<std::vector<std::string>> parameter_values = {};
    parameter_values.resize(parameter_names.size());
    if (model_runoff == "full_rainfall")
    {
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            parameter_values[0].push_back(std::to_string(cell_centre_x[idx_cell]));
            parameter_values[1].push_back(std::to_string(cell_centre_y[idx_cell]));
            parameter_values[2].push_back(std::to_string(idx_cell));
            parameter_values[3].push_back(std::to_string(precipitation_station_code[idx_cell]));
        }
    }
    else if ((model_runoff == "infiltration_coefficient") ||
             (model_runoff == "runoff_coefficient"))
    {
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            parameter_values[0].push_back(std::to_string(cell_centre_x[idx_cell]));
            parameter_values[1].push_back(std::to_string(cell_centre_y[idx_cell]));
            parameter_values[2].push_back(std::to_string(idx_cell));
            parameter_values[3].push_back(std::to_string(precipitation_station_code[idx_cell]));
            parameter_values[4].push_back(std::to_string(runoff_param [idx_cell]));
        }
    }
    else if (model_runoff == "scs")
    {
        for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
        {
            parameter_values[0].push_back(std::to_string(cell_centre_x[idx_cell]));
            parameter_values[1].push_back(std::to_string(cell_centre_y[idx_cell]));
            parameter_values[2].push_back(std::to_string(idx_cell));
            parameter_values[3].push_back(std::to_string(precipitation_station_code[idx_cell]));
            parameter_values[4].push_back(std::to_string(scs_i_a[idx_cell]));
            parameter_values[5].push_back(std::to_string(scs_s[idx_cell]));
        }
    }
    //== Create the string table for parameter values

    //-- Create file and write values
    std::string file_name = swFile::log_runoff_file;
    file_name = m_io->absolutePath(file_name);
    swIO::createFileMapString(file_name, parameter_names, parameter_values);
    //== Create file and write values
    logInfo() << "swSp - writeLogRunoff done";
}


/** @todo: adapt code */
void swDdp::writeLogSingularHeadLoss(void)
/**
 * @brief Writes log file for singular head loss
 * @par Modified attributes:
 *  - parameter_names
 *  - parameter_values
 */
{
    logInfo() << "swSp - writeSingularHeadLoss - model_singular_head_loss = " << model_singular_head_loss;
    //-- Create the string vector that will be used for the header
    std::vector<std::string> parameter_names = { "x", "y", "idx_cell" };
    parameter_names.push_back("head_loss_coeff_1");
    parameter_names.push_back("head_loss_coeff_2");
    parameter_names.push_back("head_loss_alpha");
    //== Create the string vector that will be used for the header


    //-- Create the string table for parameter values
    std::vector<std::vector<std::string>> parameter_values = {};
    parameter_values.resize(parameter_names.size());
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        parameter_values[0].push_back(std::to_string(cell_centre_x[idx_cell]));
        parameter_values[1].push_back(std::to_string(cell_centre_y[idx_cell]));
        parameter_values[2].push_back(std::to_string(idx_cell));
        parameter_values[3].push_back(std::to_string(singular_head_loss_coeff_direction_1[idx_cell]));
        parameter_values[4].push_back(std::to_string(singular_head_loss_coeff_direction_2[idx_cell]));
        parameter_values[5].push_back(std::to_string(singular_head_loss_alpha[idx_cell]));
    }
    //== Create the string table for parameter values

    //-- Create file and write values
    std::string file_name = swFile::log_singular_head_loss_file;
    file_name = m_io->absolutePath(file_name);
    swIO::createFileMapString(file_name, parameter_names, parameter_values);
    //== Create file and write values
    logInfo() << "swSp - writeSingularHeadLoss done";
}


/** @todo: adapt code */
void swDdp::writeLogSolver(void)
/**
 * @brief Writes log file for solver
 * @par No modified attribute. The solver log file is created and filled.
 */
{
    logInfo() << "swDdp - writeLogSolver";
    //-- Create file
    std::string file_name = swFile::log_solver_file;
    file_name = m_io->absolutePath(file_name);
    std::ofstream file_w;
    //== Create file

    //-- Open, Write, Close file
    file_w.open(file_name);
    //
    file_w << "cfl_max                        = " << cfl_max << "\n";
    file_w << "div_corr_option                = " << div_corr_option << "\n";
    file_w << "dt_max                         = " << dt_max << "\n";
    file_w << "dt_map                         = " << dt_map << "\n";
    file_w << "dt_time_series                 = " << dt_time_series << "\n";
    file_w << "fr_max                         = " << fr_max << "\n";
    file_w << "h_min                          = " << h_min << "\n";
    file_w << "map_storage_mode               = " << map_storage_mode << "\n";
    file_w << "model_friction                 = " << model_friction << "\n";
    file_w << "model_hyp                      = " << model_hyp << "\n";
    file_w << "model_infiltration             = " << model_infiltration << "\n";
    file_w << "model_momentum_diffusion       = " << model_momentum_diffusion << "\n";
    file_w << "model_porosity_rainfall        = " << model_porosity_rainfall << "\n";
    file_w << "model_runoff                   = " << model_runoff << "\n";
    file_w << "model_singular_head_loss       = " << model_singular_head_loss << "\n";
    file_w << "model_wind                     = " << model_wind << "\n";
    file_w << "momentum_diffusion_coefficient = " << momentum_diffusion_coefficient << "\n";
    file_w << "n_it_div_max                   = " << n_it_div_max << "\n";
    file_w << "n_por_tab_max                  = " << n_por_tab_max << "\n";
    file_w << "negative_depth_detection       = " << negative_depth_detection << "\n";
    file_w << "num_scheme                     = " << num_scheme << "\n";
    file_w << "t_0                            = " << t_0 << "\n";
    file_w << "t_max                          = " << t_max << "\n";
    file_w << "topo_discretization            = " << topo_discretization << "\n";
    file_w << "v_max                          = " << v_max << "\n";
    for (size_t idx_probe = 0; idx_probe < probe_names.size(); idx_probe++)
    {
        file_w << "probe name, x, y               = " << probe_names[idx_probe] << "\t" << 
                  probe_x[idx_probe] << "\t" << probe_y[idx_probe] << "\n";
    }

    //
    file_w.close();
    //== Open, Write, Close file

    logInfo() << "swDdp - writeLogSolver done";
}


/** @todo: adapt code */
void swDdp::writeLogWind(void)
/**
 * @brief Writes log file for wind
 * @param
 */
{
    logInfo() << "swSp - writeLogWind - model_wind = " << model_wind;
    //-- Create the string vector that will be used for the header
    std::vector<std::string> parameter_names = { "x", "y", "idx_cell", "wind_station" };
    //== Create the string vector that will be used for the header

    //-- Create the string table for parameter values
    std::vector<std::vector<std::string>> parameter_values = {};
    parameter_values.resize(parameter_names.size());
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        parameter_values[0].push_back(std::to_string(cell_centre_x[idx_cell]));
        parameter_values[1].push_back(std::to_string(cell_centre_y[idx_cell]));
        parameter_values[2].push_back(std::to_string(idx_cell));
        parameter_values[3].push_back(std::to_string(wind_station_code[idx_cell]));
    }
    //== Create the string table for parameter values

    //-- Create file and write values
    std::string file_name = swFile::log_wind_file;
    file_name = m_io->absolutePath(file_name);
    swIO::createFileMapString(file_name, parameter_names, parameter_values);
    //== Create file and write values
    logInfo() << "swSp - writeLogWind done";
}


/** @todo: adapt code */
void swDdp::writeLogs(void)
/**
 * @brief Manages calls to log writing procedures for SWES
 * @param Nothing
 * @par Modified attributes:
 *  - log files
 */
{
    //-- Write solver log
    writeLogSolver();
    //== Write solver log

    //-- Write Model logs
    writeLogBoussinesq();
    writeLogPorosity();
    writeLogEdgeFactors();

    if (model_friction != "none")
    {
        writeLogFriction();
    }
    if (model_infiltration != "none")
    {
        writeLogInfiltration();
    }
    if (model_momentum_diffusion != "none")
    {
        writeLogMomentumDiffusion();
    }
    if (model_runoff != "none")
    {
        writeLogRunoff();
    }
    if (model_singular_head_loss != "none")
    {
        writeLogSingularHeadLoss();
    }
    if (model_wind != "none")
    {
        writeLogWind();
    }
    //== Write Model logs
}


double swDdp::zOfThetaW(const size_t& idx_cell, const double& theta)
/**
 * @brief Compute the elevation z from the volume stored in the cell
 */
{
    if (theta <= 0)
    {
        return cell_z_tab[idx_cell][0];
    }
    else
    {
        size_t level = 0;
        while ((level < thetaW_tab[idx_cell].size()) && (theta > thetaW_tab[idx_cell][level]))
        {
            ++level;
        }
        double z = cell_z_tab[idx_cell][level - 1];
        z += (theta - thetaW_tab[idx_cell][level - 1]) / phiW_tab[idx_cell][level - 1]; /** @todo check division per 0 (bridge case) **/
        return z;
    }
}
