#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstring>
#include <string>
#include <vector>

#include "swSwes.h"
#include "../../../common/swFile.h"
#include "../../../common/param.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"

swSwes::swSwes(void) 
/**
 * @brief Constructor for the swSwes object
 */
{
    // Empty on purpose
}


swSwes::~swSwes(void)
/**
 * @brief Destructor for the swSwes object
 */
{
    // Empty on purpose
}

void swSwes::checkParameterConsistency(void)
/**
 * @brief checks consistency for all solver parameters
 * @return nothing, but fails if any problem is found
 */
{
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


void swSwes::computeEdgeFactors()
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
    logInfo() << "swSwes::computeEdgeFactors - Compute boundary factors";
    edge_bnd_factor_1_F1.resize(edge_bnd_nb);
    edge_bnd_factor_2_F2.resize(edge_bnd_nb);
    edge_bnd_factor_2_F3.resize(edge_bnd_nb);
    edge_bnd_factor_2_S2.resize(edge_bnd_nb);
    edge_bnd_factor_3_F2.resize(edge_bnd_nb);
    edge_bnd_factor_3_F3.resize(edge_bnd_nb);
    edge_bnd_factor_3_S2.resize(edge_bnd_nb);
    for (size_t i = 0; i < edge_bnd_nb; i++)
    {
        double w = edge_bnd_length[i];
        size_t cell_id = edge_bnd_cell_idx[i];
        double area = cell_area[cell_id];
        double nx = edge_bnd_norm_x[i];
        double ny = edge_bnd_norm_y[i];
        double coeff = w / area;

        edge_bnd_factor_1_F1[i] = + coeff;
        edge_bnd_factor_2_F2[i] = + coeff * nx;
        edge_bnd_factor_2_F3[i] = - coeff * ny;
        edge_bnd_factor_2_S2[i] = + coeff * nx;
        edge_bnd_factor_3_F2[i] = + coeff * ny;
        edge_bnd_factor_3_F3[i] = + coeff * nx;
        edge_bnd_factor_3_S2[i] = + coeff * ny;
    }
    logInfo() << "swSwes::computeEdgeFactors - Compute boundary factors done";
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
    for (size_t i = 0; i < edge_int_nb; i++)
    {
        double w = edge_int_length[i];
        size_t celL_id = edge_int_celL_idx[i];
        size_t celR_id = edge_int_celR_idx[i];
        double area_L = cell_area[celL_id];
        double area_R = cell_area[celR_id];
        double nx = edge_int_norm_x[i];
        double ny = edge_int_norm_y[i];
        double coeff_L = w / area_L;
        double coeff_R = w / area_R;

        edge_int_factor_1_F1_L[i] = - coeff_L;
        edge_int_factor_1_F1_R[i] = + coeff_R;
        edge_int_factor_2_F2_L[i] = - coeff_L * nx;
        edge_int_factor_2_F2_R[i] = + coeff_R * nx;
        edge_int_factor_2_F3_L[i] = + coeff_L * ny;
        edge_int_factor_2_F3_R[i] = - coeff_R * ny;
        edge_int_factor_2_S2_L[i] = + coeff_L * nx;
        edge_int_factor_2_S2_R[i] = + coeff_R * nx;
        edge_int_factor_3_F2_L[i] = - coeff_L * ny;
        edge_int_factor_3_F2_R[i] = + coeff_R * ny;
        edge_int_factor_3_F3_L[i] = - coeff_L * nx;
        edge_int_factor_3_F3_R[i] = + coeff_R * nx;
        edge_int_factor_3_S2_L[i] = + coeff_L * ny;
        edge_int_factor_3_S2_R[i] = + coeff_R * ny;
    }
    logInfo() << "swSwes::computeEdgeFactors - Compute internal factors done";
    //==   Internal edge factors
}


void swSwes::correctVelocity()
/**
 * @brief Eliminates abnormally large velocities
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
        double h = water_depth[idx_cell];
        if (h > h_min)
        {
            double q = unit_discharge_x[idx_cell];
            double r = unit_discharge_y[idx_cell];
            double q_norm = std::sqrt(q * q + r * r);
            double c = std::sqrt(gravity * h);
            double Fr = q_norm / c / h;
            double fact = std::max(q_norm / (h * v_max), Fr / fr_max);
            if (fact > 1.0)
            {
                limited_cell_idx.push_back(idx_cell);
                unit_discharge_x[idx_cell] = q / fact;
                unit_discharge_y[idx_cell] = r / fact;
            }
        }
    }
    //==== Eliminate abnormally large velocities
}


bool swSwes::execute()
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
    storeResultBoundaryFluxes(t);
    if (probes_nb > 0)
    {
        storeResultProbes(t);
    }
    if (has_qsection)
    {
        storeResultQsections(t);
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
        hypCompute(div_corr_option, n_it_div_max, num_scheme, dt);
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
            sourceRunoff();
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

        correctVelocity();

        ///-- Time incrementation
        t += dt;
        //== Time incrementation

        ///-- Store results if needed
        store_map = false;
        store_time_series = false;
        checkStore(t, dt);
        t_1 = clock();
        if (store_map)
        {
            std::cout << "\r  ----> calling storeResultMaps (" << t << ")                   " << std::flush;
            storeResultMaps(map_storage_mode, t);
        }
        if (store_time_series)
        {
            std::cout << "\r  ----> calling storeResultBoundaryFluxes (" << t << ")         " << std::flush;
            storeResultBoundaryFluxes(t);
        }
        if ((store_time_series) && (probes_nb > 0))
        {
            std::cout << "\r  ----> calling storeResultTimeSeries (" << t << ")             " << std::flush;
            storeResultProbes(t);
        }
        if ((store_time_series) && (has_qsection))
        {
            std::cout << "\r  ----> calling storeResultQsections (" << t << ")              " << std::flush;
            storeResultQsections(t);
        }
        while (do_extraction && (extraction_times[extraction_time_idx].first <= t))
        {
            std::cout << "\r  ----> calling storeResultExtractionTimes (" << t << ")        " << std::flush;
            storeResultExtractionTimes(map_storage_mode, t, extraction_times[extraction_time_idx].second);
        }
        total_t_clock_storage += clock() - t_1;
        //== Store results if needed

    }
    //== Time loop
    t_1 = clock();
    storeResultMaps(map_storage_mode, t);
    storeResultBoundaryFluxes(t);
    if (probes_nb > 0)
    {
        storeResultProbes(t);
        probes_output_stream.close();
    }
    if (has_qsection)
    {
        storeResultQsections(t);
        qsection_output_stream.close();
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


void swSwes::importAttributesFromSolverParam(swSolverParam *solver_param)
/**
 * @brief sets class attributes by looking into m_solver_param
 * @param[in] m_solver_param: solverParam object containing all solver parameters
 * @brief modifies required (local) computational options
 * @par Modified attributes:
 * - all swSwes.h attributes
 */
{
    biblio_reference_list = solver_param->getBiblioReferenceList(); 
    cfl_max = solver_param->getCFLMax();
    div_corr_option = solver_param->getDivergenceCorrection();
    dt_map = solver_param->getDtMap();
    dt_max = solver_param->getDtMax();
    dt_time_series = solver_param->getDtTimeSeries();
    extraction_times = solver_param->getExtractionTimes();
    fr_max = solver_param->getFrMax();
    has_qsection = solver_param->getHasQsection();
    h_min = solver_param->getHMin();
    log_write = solver_param->getLogWrite();
    map_storage_mode = solver_param->getStoreFlag();
    model_friction = solver_param->getModelFriction();
    model_hyp = solver_param->getModelHyp();
    model_infiltration = solver_param->getModelInfiltration();
    model_momentum_diffusion = solver_param->getModelMomentumDiffusion();
    model_runoff = solver_param->getModelRunoff();
    model_singular_head_loss = solver_param->getModelSingularHeadLoss();
    model_wind = solver_param->getModelWind();
    momentum_diff_max_iter = solver_param->getMomentumDiffMaxIter();
    momentum_diff_convergence_eps = solver_param->getMomentumDiffConvergenceEps();
    negative_depth_detection = solver_param->getNegativeDepthDetection();
    num_scheme = solver_param->getNumScheme();
    n_decimal = solver_param->getNDecimal();
    n_it_div_max = solver_param->getNItDivMax();
    output_format = solver_param->getOutputFormat();
    probe_names = solver_param->getProbeNames();
    probe_x = solver_param->getProbeX();
    probe_y = solver_param->getProbeY();
    topo_discretization = solver_param->getTopoDiscretization();
    t_0 = solver_param->getT0();
    t_max = solver_param->getTMax();
    v_max = solver_param->getVMax();
}


bool swSwes::init(swMesh *mesh)
/**
 * @brief Initializes the SWES model
 * @param[in] mesh: contains all mesh-related fields
 * @return bool EXIT_SUCCESS / EXIT_FAILURE
 */
{
    logInfo() << "swSwes::init - START";

    ///-- Prepare output binary file for maps
    if ((output_format == "bin") || (output_format == "all"))
    {
        std::string filename = m_io->absolutePath(swFile::output_map_binary_file);
        std::string geo_filename = m_io->absolutePath(swFile::input_geo_file);
        std::vector<std::string> param_name = {};
        if (map_storage_mode == "default")
        {
            param_name = {"water_depth", "unit_discharge_x", "unit_discharge_y", "unit_discharge_norm" };
        }
        else if (map_storage_mode == "all")
        {
            param_name = {"water_depth", "unit_discharge_x", "unit_discharge_y", "unit_discharge_norm",
                        "surface_elevation" , "velocity_x", "velocity_y", "velocity_norm" ,
                        "froude", "specific_head", "total_head" };
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
        for (size_t i = 0 ; i < probes_nb ; i++)
        {
            std::string name = probe_names[i];
            double x = probe_x[i];
            double y = probe_y[i];
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
        for (size_t i = 0 ; i < probes_nb ; i++)
        {
            std::string name = probe_names[i];
            probes_file_line += name + "_zb\t" + name + "_h\t" + name + "_q\t" + name + "_r\t";
        }
        probes_file_line += "\n";
        swIO::writeTxtLine(probes_output_stream, probes_file_line);
    }
    //==  Prepare output file for probes (first lines contain probes name and location)

    ///-- Initialize mesh   
	logInfo() << "swSwes::init - get mesh information() " ;
    setMeshTables(mesh);
	//== Initialize mesh

    ///-- Prepare output file for boundary fluxes
    std::string bnd_fluxes_filename = m_io->absolutePath(swFile::output_boundary_fluxes_file);
    bnd_fluxes_output_stream.open(bnd_fluxes_filename);
    if (bnd_fluxes_output_stream.is_open() == false)
    {
        logError() << "can't create file " << bnd_fluxes_filename;
        logFatal(__FILE__, __func__, __LINE__);
    }

    bnd_fluxes_output_stream << "time";

    for (size_t k = 0; k <list_bc_codes.size(); k++)
    {
        int code = list_bc_codes[k];
        bnd_fluxes_output_stream << "\t" << code;
    }
    bnd_fluxes_output_stream << std::endl;
    bnd_fluxes_output_stream.flush();
    //== Prepare ouput file for boundary fluxes

    ///-- Prepare output file for qsections
    std::string qsection_filename = m_io->absolutePath(swFile::output_qsection_file);
    qsection_output_stream.open(qsection_filename);
    if (qsection_output_stream.is_open() == false)
    {
        logError() << "can't create file " << qsection_filename;
        logFatal(__FILE__, __func__, __LINE__);
    }

    qsection_output_stream << "time";

    for (size_t k = 0; k <list_qsections_names.size(); k++)
    {
        std::string name = list_qsections_names[k];
        qsection_output_stream << "\t" << name;
    }
    qsection_output_stream << std::endl;
    qsection_output_stream.flush();
    //== Prepare ouput file for qsections

    ///-- Initialize variables & parameters
    ///---- Set up extraction times
    if (extraction_times.size() > 0)
    {
        do_extraction = true;
    }
    ///==== Set up extraction times

	///---- Set up initial conditions 
    // NB: variables must be declared in .h
    setInitialConditions();
    //==== Set up initial conditions

    ///---- Set up boundary conditions
    logInfo() << "Setting up boundary conditions";
    if (swParent::initBoundaryConditions() != EXIT_SUCCESS)
    {
        logError() << "Error while initialising boundary conditions";
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
    
    ///------Friction
    logInfo() << "swSwes::init - friction_model = " << model_friction;
    if (model_friction != "none")
    {
        setFrictionParameters();
    }
    //====== Friction

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
    logInfo() << "swSwes::init - compute edge factors";
    computeEdgeFactors();
    //== Compute edge factors 

    ///-- Resize fluxes and Riemann states
    h_riemann_bnd.resize(edge_bnd_nb);
    h_riemann_int.resize(edge_int_nb);
    lambda_max_bnd.resize(edge_bnd_nb);
    lambda_max_int.resize(edge_int_nb);
    lambda_minus_bnd.resize(edge_bnd_nb);
    lambda_minus_int.resize(edge_int_nb);
    lambda_plus_bnd.resize(edge_bnd_nb);
    lambda_plus_int.resize(edge_int_nb);
    q_riemann_bnd.resize(edge_bnd_nb);
    q_riemann_int.resize(edge_int_nb);
    r_riemann_bnd.resize(edge_bnd_nb);
    r_riemann_int.resize(edge_int_nb);
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
        logInfo() << "swSwes::init - write log files";
        writeLogs();
    }
    //== Write Log files (if applicable)

    logInfo() << "swSwes::init - done \n";
    return EXIT_SUCCESS;
}


double swSwes::lambdaMinus(double const& u_cell, double const& c)
/**
 * @brief Computes lambdaMinus in one cell, from u and c
 * @param[in] u_cell: velocity in the cell
 * @param[in] c: wave speed sqrt(g.h)
 * @return u_cell-c
 */
{
    return u_cell - c;
}


double swSwes::lambdaPlus(double const& u_cell, double const& c)
/**
 * @brief Computes lambdaPlus in one cell, from u and c
 * @param[in] u_cell: velocity in the cell
 * @param[in] c: wave speed sqrt(g.h)
 * @return u_cell + c
 */
{
    return u_cell + c;
}


void swSwes::setFrictionParameters(void)
/**
 * @brief read information from friction file (could be Chez, Manning, Strickler)
 * @par modified attributes:
 *  - friction_coeff_direction_1
 *  - friction_coeff_direction_2
 *  - friction_alpha
 */
{
    logInfo() << "swSwes::init - reading friction map";
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
    for (size_t i = 0; i < cell_nb; i++)
    {
        friction_coeff_direction_1[i] = param_array[0][i];
        friction_coeff_direction_2[i] = param_array[1][i];
        double alpha = param_array[2][i] * atan(1.0) / 45.0;  // Convert degrees to radians
        friction_alpha[i]     = alpha;
        friction_cos_alpha[i] = cos(alpha);
        friction_sin_alpha[i] = sin(alpha);
    }
    logInfo() << "             - friction map read \n";
}


void swSwes::setInfiltrationParameters(void)
/**
 * @brief read information from infiltration file (could be h_proportional, Horton)
 * @par modified attributes:
 *  - infiltration_param
 */
{
    logInfo() << "\n \n";
    logInfo() << "swSwes::init - reading infiltration parameter map";
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
        for (size_t i = 0; i < cell_nb; i++)
        {
            infiltration_param[i] = param_array[0][i];
        }
    }
    else if (model_infiltration == "horton")
    {
       file_name = m_io->absolutePath(swFile::input_infiltration_horton_parameter_map);
       infiltration_param.resize(cell_nb);
       swIO::readVectorParameterFile(file_name, 1, cell_nb, variable_name, param_array);
       for (size_t i = 0; i < cell_nb; i++)
       {
           infiltration_param[i] = param_array[0][i];
       }
    }
    else
    {
        logError() << "Infiltration model not handled = " << model_infiltration;
        logFatal(__FILE__, __func__, __LINE__);
    }
    logInfo() << "             - Infiltration parameter map read \n";
}


void swSwes::setInitialConditions(void)
/**
 * @brief read information from initial conditions file and assign it to state variables
 * @par modified attributes:
 *  - water_depth
 *  - unit_discharge_x
 *  - unit_discharge_y
 */
{
    logInfo() << "swSwes::init - reading initial conditions";
    std::string file_name = m_io->absolutePath(swFile::input_initial_conditions_file);

    std::vector<std::string> variable_name;
    std::vector<std::vector<double>> param_array;
    swIO::readVectorParameterFile(file_name, 3, cell_nb, variable_name, param_array);

    water_depth.resize(cell_nb);
    unit_discharge_x.resize(cell_nb);
    unit_discharge_y.resize(cell_nb);
    ///------ Assign initial water depth (if h or z present in variable names)
    bool found = false;
    auto it = std::find(variable_name.begin(), variable_name.end(), "h");
    size_t column = it - variable_name.begin();
    if (column < 3) /* Index has been found before end of vector*/
    {
        found = true;
        for (size_t i = 0; i < cell_nb; i++)
        {
            water_depth[i] = param_array[column][i];
        }
    }
    it = std::find(variable_name.begin(), variable_name.end(), "z");
    column = it - variable_name.begin();
    if (column < 3) /* Index has been found before end of vector*/
    {
        found = true;
        for (size_t i = 0; i < cell_nb; i++)
        {
            water_depth[i] = std::max(param_array[column][i] - cell_z_b[i], 0.0);
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
        for (size_t i = 0; i < cell_nb; i++)
        {
            unit_discharge_x[i] = param_array[column][i];
        }
    }
    it = std::find(variable_name.begin(), variable_name.end(), "u");
    column = it - variable_name.begin();
    if (column < 3) /* Index has been found before end of vector*/
    {
        found = true;
        for (size_t i = 0; i < cell_nb; i++)
        {
            unit_discharge_x[i] = water_depth[i] * param_array[column][i];
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
        for (size_t i = 0; i < cell_nb; i++)
        {
            unit_discharge_y[i] = param_array[column][i];
        }
    }
    it = std::find(variable_name.begin(), variable_name.end(), "v");
    column = it - variable_name.begin();
    if (column < 3) /* Index has been found before end of vector*/
    {
        found = true;
        for (size_t i = 0; i < cell_nb; i++)
        {
            unit_discharge_y[i] = water_depth[i] * param_array[column][i];
        }
    }
    if (!found)
    {
        logError() << "Could find neither r nor v in initial conditions. Stopping here.";
        logFatal(__FILE__, __func__, __LINE__);
    }
    //====== Assign initial y-unit discharge (if r or v present in variable names)
}


void swSwes::setParabParameters(void)
/**
 * @brief read information from parabolic part
 * @par modified attributes:
 *  - parab_edge_int_flux_factors 
 */
{
    /** @todo Read file and assign coefficients */
    logInfo() << "setParabParameters - top";
    //-- Prepare coefficients for matrix setup
    for (size_t idx_edge = 0; idx_edge < edge_int_nb; idx_edge++)
    {
        size_t idx_celL = edge_int_celL_idx[idx_edge];
        size_t idx_celR = edge_int_celR_idx[idx_edge];
        double dx = cell_centre_x[idx_celR] - cell_centre_x[idx_celL];
        double dy = cell_centre_y[idx_celR] - cell_centre_y[idx_celL];
        double d_norm = std::abs(dx * edge_int_norm_x[idx_edge] + dy * edge_int_norm_y[idx_edge]);
        parab_edge_int_flux_factors.push_back(edge_int_length[idx_edge] / std::max(d_norm, small_real));
    }
    //== Prepare coefficients for matrix setup
}


void swSwes::setRunoffParameters(void)
/**
 * @brief read information from runoff file (could be full rainfall, infiltration coeff, runoff coeff, SCS model)
 * @par modified attributes:
 *  - runoff_param
 *  - scs_i_a
 *  - scs_s
 */
{
    logInfo() << "\n \n";
    logInfo() << "swSwes::init - reading runoff parameter map";
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
        for (size_t i = 0; i < cell_nb; i++)
        {
            runoff_param[i] = param_array[0][i];
        }
    }
    else if (model_runoff == "runoff_coefficient")
    {
        file_name = m_io->absolutePath(swFile::input_runoff_runoff_coefficient_map);
        runoff_param.resize(cell_nb);
        swIO::readVectorParameterFile(file_name, 1, cell_nb, variable_name, param_array);
        for (size_t i = 0; i < cell_nb; i++)
        {
            runoff_param[i] = param_array[0][i];
        }
    }
    else if (model_runoff == "scs") 
    {
        //-- Read SCS model parameters
        file_name = m_io->absolutePath(swFile::input_runoff_scs_parameter_map);
        swIO::readVectorParameterFile(file_name, 2, cell_nb, variable_name, param_array);
        scs_i_a.resize(cell_nb);
        scs_s.resize(cell_nb);
        for (size_t i = 0; i < cell_nb; i++)
        {
            scs_i_a[i] = param_array[0][i];
            scs_s[i] = param_array[1][i];
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


void swSwes::setSingularHeadLossParameters(void)
/**
 * @brief read information from singular head loss file
 * @par modified attributes:
 *  - singular_head_loss_h_sing
 *  - singular_head_loss_coeff_direction_1
 *  - singular_head_loss_coeff_direction_2
 *  - singular_head_loss_coeff_alpha
 */
{
    logInfo() << "swSwes::init - reading singular head loss map";
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


void swSwes::setTransportDiffusionParameters(void)
/**
 * @brief sets transport diffusion array (cellwise)
 * @par modified attributes:
 *  - transport_scalar_diffusion_parameter
 */
{
    logInfo() << "swSwes::init - setting up transport diffusion parameters";
    transport_scalar_diffusion_parameter.resize(cell_nb);
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        transport_scalar_diffusion_parameter[idx_cell] = transport_diffusion_coefficient;
    }
    logInfo() << "             - transport diffusion parameters set up \n";
}


double swSwes::specificForce(double const& u_cell, double const& h)
/**
 * @brief swSwesHyp::specificForce
 * @param[in] u_cell: velocity in the cell
 * @param[in] h: water depth in the cell
 * @return h*u*u + 0.5*g*h*h
 */
{
    /** @todo Rewrite specific force with conserved variables? */

    return (h * u_cell * u_cell + 0.5 * gravity * h * h);
}


void swSwes::storeResultBoundaryFluxes(double const& t)
/**
 * @brief store the boundary fluxes values at the current time
 * 
 * @param[in] t current time
 */
{
    bnd_fluxes_output_stream << swString::numberToStr(t, 'e', 10);

    for (size_t k = 0; k < list_boundaries.size(); k++)
    {
        double sum = 0;
        std::vector<swEdgeBoundary *> edges = list_boundaries[k];

        for (size_t i = 0; i < edges.size(); i++)
        {
            swEdgeBoundary * edge = edges[i];
            size_t id = edge->getIdx();
            double length = edge_bnd_length[id];
            sum += F1_bnd[id] * length;
        }
        bnd_fluxes_output_stream << "\t" << swString::numberToStr(sum, 'e', 10);
    }
    bnd_fluxes_output_stream << std::endl;

    bnd_fluxes_output_stream.flush();
}


void swSwes::storeResultExtractionTimes(const std::string& map_storage_mode, double const& time, const std::string& label)
/**
 * @brief Writes the result extraction time map file
 * @param[in] map_storage_mode: string to define the output verbosity (default or all)
 * @param[in] time : current time
 * @param[in] label : label for file name
 * @par Modified attributes:
 * - content of output files
 */
{
    if (map_storage_mode != "none")
    {
        ///-- Set up file name
        std::string file_name =swFile::output_extraction_time_file_root;
        m_io->computeFileNameFromTime(t, label, file_name);
        //== Set up file name

        ///-- Fill in vectors to be written
        std::vector<std::string> param_name = {};
        if (map_storage_mode == "default")
        {
            param_name = { "x", "y", "zb", "water_depth", "unit_discharge_x", "unit_discharge_y", "unit_discharge_norm" };
        }
        else if (map_storage_mode == "all")
        {
            param_name = { "x", "y", "zb", "water_depth", "unit_discharge_x", "unit_discharge_y", "unit_discharge_norm",
                        "surface_elevation" , "velocity_x", "velocity_y", "velocity_norm" ,
                        "froude", "specific_head", "total_head" };
        }
        std::vector<std::vector<double>> param_array;
        param_array.resize(param_name.size());


        //-- Parts of the table that will always be filled in 
        for (size_t n = 0; n < cell_nb; n++) // Instructions are faster than copy
        {
            {
                param_array[0].push_back(cell_centre_x[n]);
                param_array[1].push_back(cell_centre_y[n]);
                param_array[2].push_back(cell_z_b[n]);
                param_array[3].push_back(water_depth[n]);
                param_array[4].push_back(unit_discharge_x[n]);
                param_array[5].push_back(unit_discharge_y[n]);
                double q_norm = std::sqrt(unit_discharge_x[n] * unit_discharge_x[n] + unit_discharge_y[n] * unit_discharge_y[n]);
                param_array[6].push_back(q_norm);
            }
            //== Parts of the table that will be filled in "default", "full" and "nonzero" modes
        }

        //-- Optional : only in 'all' mode
        if (map_storage_mode == "all"  )
        {
            for (size_t n = 0; n < cell_nb; n++) // Instructions are faster than copy
            {
                double h = water_depth[n];
                param_array[7].push_back(h + cell_z_b[n]);
                param_array[8].push_back(unit_discharge_x[n] / std::max(water_depth[n], small_real));
                param_array[9].push_back(unit_discharge_y[n] / std::max(water_depth[n], small_real));
                double v_norm = std::sqrt(unit_discharge_x[n] * unit_discharge_x[n] + unit_discharge_y[n] * unit_discharge_y[n])
                    / std::max(small_real, h);
                param_array[10].push_back(v_norm);
                param_array[11].push_back((v_norm) / std::max(sqrt(gravity * h), small_real));
                param_array[12].push_back(h + 0.5 * v_norm * v_norm / gravity);
                param_array[13].push_back(param_array[12][n] + cell_z_b[n]);
            }
        }
        //== Optional : only in 'all' mode

        //== Fill in vectors to be written

        ///-- Call writer
        if (swIO::createTxtFileMapDecimal(file_name, param_name, param_array, n_decimal) != EXIT_SUCCESS)
        {
            logError() << "swSwes::writeResultExtractionTimes - Exit Failure in swIO::createTxtMapDecimal \n";
            logFatal(__FILE__, __func__, __LINE__);
        }
        //== Call writer
    }

    extraction_time_idx++;
    if (extraction_time_idx >= extraction_times.size())
    {
        do_extraction = false;
    }
}


void swSwes::storeResultMaps(const std::string& map_storage_mode, double const& time)
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
            param_name = { "x", "y", "zb", "water_depth", "unit_discharge_x", "unit_discharge_y", "unit_discharge_norm" };
        }
        else if (map_storage_mode == "all")
        {
            param_name = { "x", "y", "zb", "water_depth", "unit_discharge_x", "unit_discharge_y", "unit_discharge_norm",
                        "surface_elevation" , "velocity_x", "velocity_y", "velocity_norm" ,
                        "froude", "specific_head", "total_head" };
        }
        std::vector<std::vector<double>> param_array;
        param_array.resize(param_name.size());


        //-- Parts of the table that will always be filled in 
        for (size_t n = 0; n < cell_nb; n++) // Instructions are faster than copy
        {
            {
                param_array[0].push_back(cell_centre_x[n]);
                param_array[1].push_back(cell_centre_y[n]);
                param_array[2].push_back(cell_z_b[n]);
                param_array[3].push_back(water_depth[n]);
                param_array[4].push_back(unit_discharge_x[n]);
                param_array[5].push_back(unit_discharge_y[n]);
                double q_norm = std::sqrt(unit_discharge_x[n] * unit_discharge_x[n] + unit_discharge_y[n] * unit_discharge_y[n]);
                param_array[6].push_back(q_norm);
            }
            //== Parts of the table that will be filled in "default", "full" and "nonzero" modes
        }

        //-- Optional : only in 'all' mode
        if (map_storage_mode == "all"  )
        {
            for (size_t n = 0; n < cell_nb; n++) // Instructions are faster than copy
            {
                double h = water_depth[n];
                param_array[7].push_back(h + cell_z_b[n]);
                param_array[8].push_back(unit_discharge_x[n] / std::max(water_depth[n], small_real));
                param_array[9].push_back(unit_discharge_y[n] / std::max(water_depth[n], small_real));
                double v_norm = std::sqrt(unit_discharge_x[n] * unit_discharge_x[n] + unit_discharge_y[n] * unit_discharge_y[n])
                    / std::max(small_real, h);
                param_array[10].push_back(v_norm);
                param_array[11].push_back((v_norm) / std::max(sqrt(gravity * h), small_real));
                param_array[12].push_back(h + 0.5 * v_norm * v_norm / gravity);
                param_array[13].push_back(param_array[12][n] + cell_z_b[n]);
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


void swSwes::storeResultProbes(double const& time)
/**
 * @brief Stores the probes values on one single line
 * @param[in] time : current time
 * @par Modified attributes:
 * - content of output files
 */
{   /// Line starts with time
    std::string probes_file_line = swString::numberToStr(time, 'e', 10) + "\t";
    
    ///-- For each probe, print z_b, h, q and r
    for (size_t i = 0 ; i < probes_nb ; i++)
    {
        size_t idx_cell = probe_cell_idx[i];
        probes_file_line += swString::numberToStr(cell_z_b[idx_cell], 'e', 10) + "\t" + swString::numberToStr(water_depth[idx_cell], 'e', 10) 
            + "\t" + swString::numberToStr(unit_discharge_x[idx_cell], 'e', 10) + "\t" + swString::numberToStr(unit_discharge_y[idx_cell], 'e', 10) + "\t";
    }
    probes_file_line += "\n";
    //== For each probe, print z_b, h, q and r
        
    ///-- Call writer to print line in probe file
    if (swIO::writeTxtLine(probes_output_stream, probes_file_line) != EXIT_SUCCESS)
    {
        logError() << "swSwes::storeResultProbes - Exit Failure in swIO::writeTxtLine \n";
        logFatal(__FILE__, __func__, __LINE__);
    }
    //== Call writer
}


void swSwes::storeResultQsections(double const& t)
/**
 * @brief store the qsections values at the current time
 * 
 * @param[in] t current time
 */
{
    qsection_output_stream << swString::numberToStr(t, 'e', 10);

    for (size_t k = 0; k < list_qsections.size(); k++)
    {
        double sum = 0;
        std::vector<std::string> edges = list_qsections[k];
        std::vector<int> signs = list_qsections_signs[k];

        for (size_t i = 0; i < edges.size(); i++)
        {
            std::vector<std::string> edge = swString::splitStr(edges[i],"_");
            std::string type = edge[0];
            bool conversion_success = true;
            size_t id = swString::stringToUInt(edge[1], conversion_success);
            if (type == "bnd")
            {
                double length = edge_bnd_length[id];
                sum += F1_bnd[id] * length * signs[i];
            }
            else if (type == "int")
            {
                double length = edge_int_length[id];
                sum += F1_int[id] * length * signs[i];
            }
        }
        qsection_output_stream << "\t" << swString::numberToStr(sum, 'e', 10);
    }
    qsection_output_stream << std::endl;

    qsection_output_stream.flush();
}


void swSwes::writeLogEdgeFactors(void)
/**
 * @brief Writes log file for edge factors
 * @par Modified attributes:
 */
{
    logInfo() << "swSwes - writeLogEdgeFactors";

    //-- Boundary edges
    //---- Create the string vector that will be used for the header
    std::vector<std::string> parameter_names = {"x", "y", "idx_edge"};
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
    for (size_t idx_edge = 0; idx_edge < edge_bnd_nb; idx_edge++)
    {
        parameter_values[0].push_back(std::to_string(edge_bnd_midpoint_x[idx_edge]));
        parameter_values[1].push_back(std::to_string(edge_bnd_midpoint_y[idx_edge]));
        parameter_values[2].push_back(std::to_string(idx_edge));
        parameter_values[3].push_back(std::to_string(edge_bnd_factor_1_F1[idx_edge]));
        parameter_values[4].push_back(std::to_string(edge_bnd_factor_2_F2[idx_edge]));
        parameter_values[5].push_back(std::to_string(edge_bnd_factor_2_F3[idx_edge]));
        parameter_values[6].push_back(std::to_string(edge_bnd_factor_2_S2[idx_edge]));
        parameter_values[7].push_back(std::to_string(edge_bnd_factor_3_F2[idx_edge]));
        parameter_values[8].push_back(std::to_string(edge_bnd_factor_3_F3[idx_edge]));
        parameter_values[9].push_back(std::to_string(edge_bnd_factor_3_S2[idx_edge]));
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
    parameter_names = {"x", "y", "idx_edge"};
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
    for (size_t idx_edge = 0; idx_edge < edge_int_nb; idx_edge++)
    {
        parameter_values[0].push_back(std::to_string(edge_int_midpoint_x[idx_edge]));
        parameter_values[1].push_back(std::to_string(edge_int_midpoint_y[idx_edge]));
        parameter_values[2].push_back(std::to_string(idx_edge));
        parameter_values[3].push_back(std::to_string(edge_int_factor_1_F1_L[idx_edge]));
        parameter_values[4].push_back(std::to_string(edge_int_factor_1_F1_R[idx_edge]));
        parameter_values[5].push_back(std::to_string(edge_int_factor_2_F2_L[idx_edge]));
        parameter_values[6].push_back(std::to_string(edge_int_factor_2_F2_R[idx_edge]));
        parameter_values[7].push_back(std::to_string(edge_int_factor_2_F3_L[idx_edge]));
        parameter_values[8].push_back(std::to_string(edge_int_factor_2_F3_R[idx_edge]));
        parameter_values[9].push_back(std::to_string(edge_int_factor_2_S2_L[idx_edge]));
        parameter_values[10].push_back(std::to_string(edge_int_factor_2_S2_R[idx_edge]));
        parameter_values[11].push_back(std::to_string(edge_int_factor_3_F2_L[idx_edge]));
        parameter_values[12].push_back(std::to_string(edge_int_factor_3_F2_R[idx_edge]));
        parameter_values[13].push_back(std::to_string(edge_int_factor_3_F3_L[idx_edge]));
        parameter_values[14].push_back(std::to_string(edge_int_factor_3_F3_R[idx_edge]));
        parameter_values[15].push_back(std::to_string(edge_int_factor_3_S2_L[idx_edge]));
        parameter_values[16].push_back(std::to_string(edge_int_factor_3_S2_R[idx_edge]));
    }
    //==== Create the string table for parameter values

    //---- Create file and write values
    file_name = swFile::log_edge_bnd_factors_file;
    file_name = m_io->absolutePath(file_name);
    swIO::createFileMapString(file_name, parameter_names, parameter_values);
    //==== Create file and write values
    //== Boundary edges


    logInfo() << "swSwes - writeLogEdgeFactors done";
}


void swSwes::writeLogFriction(void)
/**
 * @brief Writes log file for friction
 * @par Modified attributes:
 *  - parameter_names
 *  - parameter_values
 */
{
    logInfo() << "swSwes - writeLogFriction - model_friction = " << model_friction;
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
    logInfo() << "swSwes - writeLogFriction done";
}


void swSwes::writeLogInfiltration(void)
/**
 * @brief Writes log file for infiltration
 * @par Modified attributes:
 *  - parameter_names
 *  - parameter_values
 */
{
    logInfo() << "swSwes - writeLogInfiltration - model_infiltration = " << model_infiltration;
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
    logInfo() << "swSwes - writeLogInfiltration done";
}


void swSwes::writeLogMomentumDiffusion(void)
/**
 * @brief Writes log file for runoff
 * @param
 * @todo Write the code
 */
{
    logInfo() << "Not implemented yet";
}


void swSwes::writeLogRunoff(void)
/**
 * @brief Writes log file for runoff
 * @par Modified attributes:
 *  - parameter_names
 *  - parameter_values
 */
{
    logInfo() << "swSwes - writeLogRunoff - model_runoff = " << model_runoff;
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
    logInfo() << "swSwes - writeLogRunoff done";
}


void swSwes::writeLogSingularHeadLoss(void)
/**
 * @brief Writes log file for singular head loss
 * @par Modified attributes:
 *  - parameter_names
 *  - parameter_values
 */
{
    logInfo() << "swSwes - writeSingularHeadLoss - model_singular_head_loss = " << model_singular_head_loss;
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
    logInfo() << "swSwes - writeSingularHeadLoss done";
}


void swSwes::writeLogSolver(void)
/**
 * @brief Writes log file for solver
 * @par No modified attribute. The solver log file is created and filled.
 */
{
    logInfo() << "swSwes - writeLogSolver";
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
    file_w << "model_runoff                   = " << model_runoff << "\n";
    file_w << "model_singular_head_loss       = " << model_singular_head_loss << "\n";
    file_w << "model_wind                     = " << model_wind << "\n";
    file_w << "momentum_diffusion_coefficient = " << momentum_diffusion_coefficient << "\n";
    file_w << "n_it_div_max                   = " << n_it_div_max << "\n";
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

    logInfo() << "swSwes - writeLogSolver done";
}


void swSwes::writeLogWind(void)
/**
 * @brief Writes log file for wind
 * @param
 */
{
    logInfo() << "swSwes - writeLogWind - model_wind = " << model_wind;
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
    logInfo() << "swSwes - writeLogWind done";
}


void swSwes::writeLogs(void)
/**
 * @brief Manages calls to log writing procedures for SWES
 * @param Nothing
 */
{
    //-- Write solver log
    writeLogSolver();
    //== Write solver log

    //-- Write Model logs
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
