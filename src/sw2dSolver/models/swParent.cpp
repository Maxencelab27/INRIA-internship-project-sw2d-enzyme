#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <unordered_map>

#include "swParent.h"
#include "./ddp/swDdp.h"
#include "./dip/swDip.h"
#include "./liswes/swLiswes.h"
#include "./sp/swSp.h"
#include "./swes/swSwes.h"

#include "../../common/swFile.h"
#include "../../common/geometry/swMesh.h"
#include "../../common/geometry/swCell.h"
#include "../../common/geometry/swEdgeBoundary.h"
#include "../../common/geometry/swEdgeIntern.h"
#include "../../common/geometry/swNode.h"
#include "../../common/swIO.h"
#include "../../common/swLog.h"
#include "../../common/param.h"
#include "common/swBoundaryConditionRatingCurve.h"
#include "common/swBoundaryConditionTimeSeries.h"
#include "common/swSource.h"

swParent::swParent(void)
/**
 * @brief Construct the swParent object
 *
 */
{
    // Empty
}


swParent::~swParent(void)
/**
 * @brief Destroy the swParent object
 *
 */
{
    //m_mesh = nullptr;
}


void swParent::checkStore(const double& time, const double& delta_time)
/**
 * @brief sts the value of parameters store_map and store_time_series
 * @param[in] time (double) current computational time (at the end of time step)
 * @param[in] delta_time (double) time step
 * @par Modified attributes :
 *  - store_map
 *  - store_time_series
 * @return Nothing
 */
{
    //-- Map storage-related time step limitation
    if (dt_map <= -small_real) 
    {
        // Negative dt means maps will be stored every time step
        store_map = true;
    }
    else if (dt_map > small_real)
    {
        if ((time/dt_map) == (int)(time/dt_map)) /* Storing is mandatory */
        {
            store_map = true;
        }
    }
    //== Map storage-related time step limitation

    //-- Determine whether time series should be stored
    if (dt_time_series <= -small_real)
    {
        // Negative dt_time_series time series will be stored every time step
        store_time_series = true;
    }
    else if (dt_time_series > small_real)
    {
        int t_now_over_dt_int = (int)((time - delta_time) / dt_time_series);
        int t_later_over_dt_int = (int)(time / dt_time_series);
        if (t_later_over_dt_int > t_now_over_dt_int) /* Storing is mandatory */
        {
            store_time_series = true;
        }
    }
    //== Determine whether time series should be stored
}


void swParent::computeBC(double const& t)
/**
 * @brief Fills in the boundary bc_val_interpolated vector with numerical values
 * @param[in] t Current simulation timeNothing
 * @par Modified attributes:
 *   - bc_val_interpolated
 */
{
    /** @todo Implement for real */
    ///-- Loop over all BC types
    for (size_t bc_idx = 0; bc_idx < bc_type_nb; bc_idx++)
    {
        if (bc_type[bc_idx] == "rc1")
        // Power law-type rating curve
        {
            logWarn() << "RC1 NOT IMPLEMENTED YET !!!";
        }
        else if (bc_type[bc_idx] == "rc2")
            // Interpolated rating curve
        {
            logWarn() << "RC2 NOT IMPLEMENTED YET !!!";
        }
        else
        // BC is provided in the form of time series that must be interpolated from the file
        /// @todo Remove access to ->time_values and ->bc_values within the loop 
        {
            bc_val_interpolated.push_back(interpolateTimeSeries(bc_timeseries_list[bc_idx]->time_values, bc_timeseries_list[bc_idx]->bc_values, t));
        }
    }
    //==  Loop over all BC types
}



swParent* swParent::create(swSolverParam *solver_param, swIO *io)
/**
 * @brief Create the swParent object
 * @param[in] solver_param : object containg all solver parameters
 * @return swParent object
 */
{
    swParent *model = nullptr;
    logInfo() << "  model_hyp = " << solver_param->getModelHyp();
    if (solver_param->getModelHyp() == "swes")
    {
        model = new swSwes();
    }
    else if (solver_param->getModelHyp() == "ddp")
    {
        model = new swDdp();
    }
    else if (solver_param->getModelHyp() == "dip")
    {
        model = new swDip();
    }
    else if (solver_param->getModelHyp() == "liswes")
    {
        model = new swLiswes();
    }
    /*
    else if (solver_param->getModelHyp() == "mp")
    {
        model = new swMp();
    }
    */
    else if (solver_param->getModelHyp() == "sp")
    {
        model = new swSp();
    }
    else
    {
        logError() << "Hyperbolic model not implemented:" << solver_param->getModelHyp();
        logError() << "  Allowed hyperbolic models (keyword model_hyperb) are:";
        logError() << "    - dip    : dual Integral Porosity model";
        logError() << "    - liswes : local inertial shallow water model";
        //logError() << "    - mp     : multiple Porosity model";
        logError() << "    - sp     : single porosity model";
        logError() << "    - swes   : classical shallow water model";
        logFatal(__FILE__, __func__, __LINE__);
    }
    model->setIOPointer(io);    // pass the pointer to the swIO object to the model
    return model;
}


void swParent::displayBiblio(void)
/**
 * @brief Diplays the references at the end of simulation
 * @param none
 */
{
    logInfo();
    logInfo() << "*** REFERENCES ***";
    logInfo() << "The SW2D modules you just used are based on (and not limited to) the research publications hereafter";
    for (size_t i = 0; i < biblio_reference_list.size(); i++)
    {
        logInfo() << i + 1 << ") " << biblio_reference_list[i];
    }

}

bool swParent::initBoundaryConditions()
/**
 * @brief initialize boundary conditions (read input file)
 * @return bool EXIT_SUCCESS / EXIT_FAILURE;
 */
{
    std::string filename;

    logInfo() << "  Assigning the boundary time series";
    filename = m_io->absolutePath(swFile::input_hydro_boundary_timeseries);
    std::vector<int> bc_timeseries_codes;
    std::vector<std::string> bc_timeseries_type;
    std::vector<double> timebc;
    std::vector<std::vector<double>> val_bc;
    if (m_io->readBoundaryTimeSeriesFile(filename,
                                                  bc_timeseries_codes,
                                                  bc_timeseries_type,
                                                  timebc, val_bc) != EXIT_SUCCESS)
    {
        logError() << "Error while reading boundary time series";
        return EXIT_FAILURE;
    }

    for (size_t i = 0; i < bc_timeseries_codes.size(); ++i)
    {
        swBoundaryConditionTimeSeries *bc = new swBoundaryConditionTimeSeries(bc_timeseries_type[i]);
        bc->initialise(timebc, val_bc[i]);
        bc_timeseries_list.push_back(bc);
    }
    logInfo() << "  ... Done";


    logInfo() << "  Assigning the rating curves";
    filename = m_io->absolutePath(swFile::input_hydro_boundary_rating_curve);
    std::vector<int> bc_ratingcurve_codes;
    if (std::filesystem::exists(filename) == true)
    {
        //-- read the BC rating curve data
        std::vector<std::vector<std::string>> rating_curve_data = {};
        if (m_io->readBoundaryRatingCurveFile(filename, rating_curve_data) != EXIT_SUCCESS)
        {
            logError() << "Error while reading boundary rating curves";
            return EXIT_FAILURE;
        }
        //-- read the BC rating curve data

        for (size_t i = 0; i < rating_curve_data.size(); ++i)
        {
            std::vector<std::string> lines = rating_curve_data[i];
            // lines[0] -> BC_code  -4
            // lines[1] -> Type  interp
            // lines[2] -> z0   0.3
            // lines[3,4] -> dépend

            //-- getting bc type and creating the corresponding bc
            std::string line = lines[1];
            std::vector<std::string> list = swString::splitStrOnce(swString::trim(line), swString::regexWhiteSpace(), 1);
            if ((list.size() != 2) || (list[0].compare("type") != 0))
            {
                logError() << "Impossible to create rating curve boundary - 1st parameter";
                logError() << lines;
                logFatal(__FILE__, __func__, __LINE__);
            }
            std::string bc_type = list[1];

            swBoundaryConditionRating *bc = swBoundaryConditionRating::create(bc_type);
            //== getting bc type and creating the corresponding bc

            //-- initialise the rating bc
            bc->initialise(lines);
            bc_ratingcurve_list.push_back(bc);
            //== initialise the rating bc

            //-- storing the bc code for the rating curve
            line = lines[0];
            list = swString::splitStrOnce(swString::trim(line), swString::regexWhiteSpace(), 1);
            if ((list.size() != 2) || (list[0].compare("bc_code") != 0))
            {
                logError() << "Impossible to create rating curve boundary - 1st parameter";
                logError() << lines;
                logFatal(__FILE__, __func__, __LINE__);
            }
            bool conversion_success = true;
            int bc_code = swString::stringToDouble(list[1], conversion_success);
            if ((conversion_success == false) || (bc_code >= 0))
            {
                logError() << "Impossible to create rating curve boundary - 1st parameter";
                logError() << "  Parameter 'bc_code' has to be a negative integer";
                logError() << lines;
                logFatal(__FILE__, __func__, __LINE__);
            }
            bc_ratingcurve_codes.push_back(bc_code);
            //== storing the bc code for the rating curve
        }
        logInfo() << "  ... Done";
    }
    else
    {
        logInfo() << "  ... no rating curve boundary condition";
    }

    std::unordered_map<int, swBoundaryCondition*> bc_map;   /**< bc_map is a table containing:
    the bc_code (an int) and a pointer to the corresponding bc
    It is used to provide the pointer to the correct bc class object to the swEdgeBoundary object
    based on the bc_code (stored attributes of the swEdgeBoundary class).
    The use of this variable is local to this function
    **/
    for (size_t k = 0; k < bc_timeseries_codes.size(); ++k)
    {
        bc_map.insert({bc_timeseries_codes[k], bc_timeseries_list[k]});
    }
    for (size_t k = 0; k < bc_ratingcurve_codes.size(); ++k)
    {
        bc_map.insert({bc_ratingcurve_codes[k], bc_ratingcurve_list[k]});
    }

    for (size_t i = 0; i < edge_bnd_nb; ++i)
    {
        int code = edge_bnd_bc_code[i];
        swBoundaryCondition *bc = nullptr;
        if (bc_map.find(code) == bc_map.end())
        {
            logError() << "Error while assigning the boundary condtions";
            logError() << "Looking for bc code:" << code;
            logError() << " not found in " << swFile::input_hydro_boundary_rating_curve
                       << "or" << swFile::input_hydro_boundary_timeseries;
            logFatal(__FILE__, __func__, __LINE__);
        }
        else
        {
            bc = bc_map[code];
        }

        edge_bnd_bc_pointer.push_back(bc);
    }
    return EXIT_SUCCESS;
}


double swParent::interpolateTimeSeries(std::vector<double> const& times_vector, std::vector<double> const& vals_vector,
                                       double const& t) 
/**
 * @brief Interpolates linearly within time series from a table containing numerical values
 * @param[in] times_vector[n] : time for the n-th record
 * @param[in] vals_vector[n]  : variable value for the n-th record
 * @param[in] t : the time at which the variable is to be interpolated
 * @return Interpolated value
*/
{
    size_t nb_times = times_vector.size();
    size_t nb_vals  = vals_vector.size();

    //-- Check sizes
    if (nb_times != nb_vals)
    {
        logError() << "Time vector and Values vector have different sizes - stopping";
        logError() << "Time vector size   = " << nb_times;
        logError() << "Values vector size = " << nb_vals;
        logFatal(__FILE__, __func__, __LINE__);
    }
    //== Check sizes

    //-- Configurations where no interpolation is needed
    if (t <= times_vector[0] - small_real)
    {
        return vals_vector[0];
    }
    else if (t >= times_vector[nb_times - 1])
    {
        return vals_vector[nb_times - 1];
    }
    //== Configurations where no interpolation is needed
    
    //-- Otherwise: interpolate
    else 
    {
        size_t iterator = 0;
        double t_1 = times_vector[iterator];
        double t_2 = times_vector[iterator + 1];
        double v_1 = vals_vector[iterator];
        double v_2 = vals_vector[iterator + 1];

        while ((t > t_2) && (iterator < nb_times - 1))
        // sweep time vector until $t \in [t_1 , t_2]$
        {
            v_1 = v_2;
            t_1 = t_2;
            iterator++;
            t_2 = times_vector[iterator];
            v_2 = vals_vector[iterator];
        }
        // Linear interpolation
        double val = v_1 + (v_2 - v_1) * (t - t_1) / (t_2 - t_1);
        return val;
    }
    //== Otherwise: interpolate
}


void swParent::setMeshTables(swMesh *mesh)
/**
 * @brief Casts the attributes of sMesh*, swCell*, swEdge* classes to tables for the sake of CPU efficiency
 * @param[in] mesh : mesh class objerct containing all geometry information
 * @par Modified attributes:
 *   - all SWES attributes related to heometry
 */
{
    logInfo() << "  swParent::setMeshTables - top " ;
    cell_nb = mesh->getCellsList().size();   /**< number of cells **/
    edge_bnd_nb = mesh->getEdgeBoundaryList().size();   /**< number of boundary edges **/
    edge_int_nb = mesh->getEdgeInternList().size();   /**< number of internal edges **/
    node_nb = mesh->getNodesList().size();       /**< number of nodes **/
    logInfo() << "  swParent::setMeshTables - number of cells, nodes, edges_bnd, edges_int = "
	          << cell_nb << " "	<< node_nb << " " << edge_bnd_nb << " " << edge_int_nb ;
	
    //-- Fill in node-based tables
    for (swNode* node : mesh->getNodesList())
	{
		node_idx.push_back(node->getIdx());
		node_coord_x.push_back(node->getX());
		node_coord_y.push_back(node->getY());
		node_coord_z.push_back(node->getZ());
	}
    logInfo() << "  nodes ...done";
	//== Fill in node-based tables
	
	//-- Fill in boundary edge-based tables
    for (swEdgeBoundary* edge_bnd : mesh->getEdgeBoundaryList())
    {
        edge_bnd_cell_idx.push_back(edge_bnd->getCell()->getIdx());
        edge_bnd_idx.push_back(edge_bnd->getIdx());
        edge_bnd_length.push_back(edge_bnd->getLength());
        edge_bnd_midpoint_x.push_back(edge_bnd->getXM());
        edge_bnd_midpoint_y.push_back(edge_bnd->getYM());
        edge_bnd_midpoint_z.push_back(edge_bnd->getZM());
        edge_bnd_node_begin_x.push_back(edge_bnd->getNodes()[0]->getX());
        edge_bnd_node_begin_y.push_back(edge_bnd->getNodes()[0]->getY());
        edge_bnd_node_end_x.push_back(edge_bnd->getNodes()[1]->getX());
        edge_bnd_node_end_y.push_back(edge_bnd->getNodes()[1]->getY());
        edge_bnd_norm_x.push_back(edge_bnd->getNormX());
        edge_bnd_norm_y.push_back(edge_bnd->getNormY());
        edge_bnd_bc_code.push_back(edge_bnd->getBCCode());
    }
    logInfo() << "  boundary edges ...done";
    //== Fill in boundary edge-based tables

    //-- Fill in internal edge-based tables
    for (swEdgeIntern* edge_int : mesh->getEdgeInternList())
    {
        edge_int_celL_idx.push_back(edge_int->getCells()[0]->getIdx());
        edge_int_celR_idx.push_back(edge_int->getCells()[1]->getIdx());
        edge_int_idx.push_back(edge_int->getIdx());
        edge_int_length.push_back(edge_int->getLength());
        edge_int_midpoint_x.push_back(edge_int->getXM());
        edge_int_midpoint_y.push_back(edge_int->getYM());
        edge_int_midpoint_z.push_back(edge_int->getZM());
        edge_int_node_begin_x.push_back(edge_int->getNodes()[0]->getX());
        edge_int_node_begin_y.push_back(edge_int->getNodes()[0]->getY());
        edge_int_node_end_x.push_back(edge_int->getNodes()[1]->getX());
        edge_int_node_end_y.push_back(edge_int->getNodes()[1]->getY());
        edge_int_norm_x.push_back(edge_int->getNormX());
        edge_int_norm_y.push_back(edge_int->getNormY());
    }
    logInfo() << "  internal edges ...done";
    //== Fill in boundary edge-based tables

    //-- Fill in cell-based tables
    for (swCell* cell : mesh->getCellsList())
    {
        cell_area.push_back(cell->getArea());
        cell_centre_x.push_back(cell->getXG());
        cell_centre_y.push_back(cell->getYG());
        cell_id.push_back(cell->getIdx());
        cell_perimeter.push_back(cell->getPerimeter());
        cell_z_b.push_back(cell->getZBottom());
    }

    cell_edges_bnd_idx.resize(cell_nb);
    for (size_t i = 0; i < edge_bnd_nb; i++)
    {
        size_t cell_idx = edge_bnd_cell_idx[i];
        cell_edges_bnd_idx[cell_idx].push_back(i);
    }

    cell_edges_int_idx.resize(cell_nb);
    cell_neigh_cell_idx.resize(cell_nb);
    for (size_t i = 0; i < edge_int_nb; i++)
    {
        size_t celL_idx = edge_int_celL_idx[i];
        size_t celR_idx = edge_int_celR_idx[i];
        cell_edges_int_idx[celL_idx].push_back(i);
        cell_edges_int_idx[celR_idx].push_back(i);
        cell_neigh_cell_idx[celL_idx].push_back(celR_idx);
        cell_neigh_cell_idx[celR_idx].push_back(celL_idx);
    }
    logInfo() << "  cells ...done";
    //== Fill in cell-based tables

    //-- Fill in q sections tables
    list_qsections = mesh->getQSections();
    list_qsections_names = mesh->getQSectionNames();
    list_qsections_signs = mesh->getQSectionSigns(); 
    //== Fill in q sections tables

    //-- Fill in boundaries tables
    list_bc_codes = mesh->getBcCodes();
    list_boundaries = mesh->getBoundaries();
    //== Fill in boundaries tables
    
}


void swParent::setIOPointer(swIO *io)
/**
 * @brief set internal attribute of type swIO
 * @param[in] io : pointer to a swIO object
 * @par Modified attributes :
 *   - m_io
 */
{
    m_io = io;
}


void swParent::setPrecipitationStationCodes(void)
/**
 * @brief read precipitation station codes values 
 * @par Modified attributes:
 *   - precipitation_station_code
 */
{

    std::string file_name = m_io->absolutePath(swFile::input_precipitation_stations);
    std::vector<std::string> variable_name;
    std::vector<std::vector<double>> param_array;
    swIO::readVectorParameterFile(file_name, 1, cell_nb, variable_name, param_array);
    precipitation_station_code.resize(cell_nb);
    double code_integer_part;
    double code_fractional_part;
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        code_fractional_part = std::abs(std::modf(param_array[0][idx_cell], &code_integer_part));
        if ((code_fractional_part > small_real) || (param_array[0][idx_cell] < 0))
        {
            logError() << "Non-integer/negative wind code read for cell " << idx_cell << " - value = " << param_array[0][idx_cell];
            logFatal(__FILE__, __func__, __LINE__);
        }
        precipitation_station_code[idx_cell] = param_array[0][idx_cell];
    }
}


void swParent::setPrecipitationTimeSeriesTable(void)
/**
 * @brief Sets the tables from which the precipitation time series will be interpolated
 * @param Nothing
 * @par Modified attributes:
 *  - precipitation_times_vector
 *  - precipitation_rate_table
 */
{
    std::string file_name = m_io->absolutePath(swFile::input_precipitation_time_series);
    swIO::readTimeSeriesScalarFile(file_name, precipitation_times_vector, stationwise_precipitation_rate_table);

    for (size_t n = 0; n < precipitation_times_vector.size(); n++)
    {
        logInfo() << "    " << precipitation_times_vector[n];
    }
    for (size_t idx_station = 0; idx_station < stationwise_precipitation_rate_table.size(); idx_station++)
    {
        logInfo() << "  Station " << idx_station;
        for (size_t n = 0; n < wind_times_vector.size(); n++)
        {
            logInfo() << "    precip. rate = " << stationwise_precipitation_rate_table[idx_station][n];
        }
    }
}


void swParent::setWindStationCodes(void)
/**
 * @brief read wind station codes values 
 * @par Modified attributes:
 *   - wind_station_code
 */
{
    std::string file_name = m_io->absolutePath(swFile::input_wind_stations);
    std::vector<std::string> variable_name;
    std::vector<std::vector<double>> param_array;
    swIO::readVectorParameterFile(file_name, 1, cell_nb, variable_name, param_array);
    wind_station_code.resize(cell_nb);
    double code_integer_part;
    double code_fractional_part;
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        code_fractional_part = std::abs(std::modf(param_array[0][idx_cell], &code_integer_part));
        if ((code_fractional_part > small_real) || (param_array[0][idx_cell] < 0))
        {
            logError() << "Non-integer/negative wind code read for cell " << idx_cell << " - value = " << param_array[0][idx_cell];
            logFatal(__FILE__, __func__, __LINE__);
        }
        wind_station_code[idx_cell] = param_array[0][idx_cell];
    }
}


void swParent::setWindVelocityTimeSeriesTables(void)
/**
 * @brief Sets the tables from which the wind time series will be interpolated
 * @param Nothing
 * @par Modified attributes:
 *  - wind_velocity_x_table
 *  - wind_velocity_y_table
 */
{
    std::string file_name = m_io->absolutePath(swFile::input_wind_time_series);
    std::vector<std::vector<double>> norm_array = {};
    std::vector<std::vector<double>> angle_degrees_array = {};
    swIO::readTimeSeries2dVectorFile(file_name, wind_times_vector, norm_array, angle_degrees_array);
    wind_velocity_x_table.resize(norm_array.size());
    wind_velocity_y_table.resize(norm_array.size());
    logInfo() << "setWindVelocityTimeSeries() - calling swSource::windConvertFromPolarToCartesian";
    for (size_t idx_station = 0; idx_station < norm_array.size(); idx_station++)
    {
        swSource::windConvertFromPolarToCartesian(norm_array[idx_station], angle_degrees_array[idx_station],
            wind_velocity_x_table[idx_station], wind_velocity_y_table[idx_station]);
    }
}
