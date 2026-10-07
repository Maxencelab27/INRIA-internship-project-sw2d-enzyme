#include "swIO.h"

#include <filesystem>
#include <regex>
#include <fstream>

#include "swFile.h"
#include "../common/swLog.h"
#include "../common/swString.h"


swIO::swIO(const std::string &user_defined_folder)
/**
 * @brief Constructor 
 * @param[in] user_defined_folder (string) the folder defined by the user
 * @par Modified attribute: project_folder
 */
{
    const std::string last_char = &user_defined_folder[user_defined_folder.size() - 1];
    if (last_char.compare(swString::dirSeparator()) == 0)
    {
        // user_defined_folder ended by directory separator
        project_folder = user_defined_folder;
    }
    else
    {
        // user_defined_folder not ended by directory separator
        project_folder = user_defined_folder + swString::dirSeparator();
    }
}

swIO::~swIO(void)
/**
 * @brief Destructor 
 */
{
}


std::string swIO::absolutePath(const std::string &file)
/**
 * @brief swIO::absolutePath provides the complete path of a file (relative to where the executable was called)
 * @param[in] file
 * @return the path of a file/folder (relative to where the executable was called)
 */
{
    return project_folder + file;
}


void swIO::computeFileNameFromTime(double const& t_time, std::string& file_name)
/**
 * @brief Sets up the name of the Hydro result file from the simulation time in seconds
 * @param[in] t_time Simulation time
 * @param[in, out] file_name of the file
*/
{
    size_t nb_d = int(t_time / 86400.);
    size_t nb_d_e3 = int(nb_d / 1000.);
    size_t nb_d_e2 = int((nb_d - 1000. * nb_d_e3) / 100.);
    size_t nb_d_e1 = int((nb_d - 1000. * nb_d_e3 - 100. * nb_d_e2) / 10.);
    size_t nb_d_e0 = int(nb_d - 1000. * nb_d_e3 - 100. * nb_d_e2 - 10. * nb_d_e1);
    size_t nb_h = int((t_time - float(nb_d) * 86400.) / 3600.);
    size_t nb_h_e1 = int(nb_h / 10.);
    size_t nb_h_e0 = int(nb_h - 10. * nb_h_e1);
    size_t nb_m = int((t_time - nb_d * 86400. - nb_h * 3600.) / 60.);
    size_t nb_m_e1 = int(nb_m / 10.);
    size_t nb_m_e0 = int(nb_m - 10. * nb_m_e1);
    size_t nb_s = int(t_time - 86400. * nb_d - 3600. * nb_h - 60. * nb_m);
    size_t nb_s_e1 = int(nb_s / 10.);
    size_t nb_s_e0 = int(nb_s - 10. * nb_s_e1);
    size_t nb_mill = int((t_time - int(t_time)) * 1000. + 0.499999);
    size_t nb_mill_e2 = int(nb_mill / 100.);
    size_t nb_mill_e1 = int((nb_mill - 100 * nb_mill_e2) / 10.);
    size_t nb_mill_e0 = int(nb_mill - 100 * nb_mill_e2 - 10 * nb_mill_e1);

    file_name += std::to_string(nb_d_e3)
              + std::to_string(nb_d_e2) + std::to_string(nb_d_e1)
              + std::to_string(nb_d_e0) + "d"
              + std::to_string(nb_h_e1) + std::to_string(nb_h_e0) + "h"
              + std::to_string(nb_m_e1) + std::to_string(nb_m_e0) + "m"
              + std::to_string(nb_s_e1) + std::to_string(nb_s_e0) + "s"
              + std::to_string(nb_mill_e2) + std::to_string(nb_mill_e1)
              + std::to_string(nb_mill_e0) + ".txt";
    file_name = absolutePath(file_name);
}


void swIO::computeFileNameFromTime(double const& t_time, std::string label, std::string& file_name)
/**
 * @brief Sets up the name of the Hydro result file from the simulation time in seconds
 * @param[in] t_time Simulation time
 * @param[in] label Label for extraction time files
 * @param[in, out] name of the file
*/
{
    size_t nb_d = int(t_time / 86400.);
    size_t nb_d_e3 = int(nb_d / 1000.);
    size_t nb_d_e2 = int((nb_d - 1000. * nb_d_e3) / 100.);
    size_t nb_d_e1 = int((nb_d - 1000. * nb_d_e3 - 100. * nb_d_e2) / 10.);
    size_t nb_d_e0 = int(nb_d - 1000. * nb_d_e3 - 100. * nb_d_e2 - 10. * nb_d_e1);
    size_t nb_h = int((t_time - float(nb_d) * 86400.) / 3600.);
    size_t nb_h_e1 = int(nb_h / 10.);
    size_t nb_h_e0 = int(nb_h - 10. * nb_h_e1);
    size_t nb_m = int((t_time - nb_d * 86400. - nb_h * 3600.) / 60.);
    size_t nb_m_e1 = int(nb_m / 10.);
    size_t nb_m_e0 = int(nb_m - 10. * nb_m_e1);
    size_t nb_s = int(t_time - 86400. * nb_d - 3600. * nb_h - 60. * nb_m);
    size_t nb_s_e1 = int(nb_s / 10.);
    size_t nb_s_e0 = int(nb_s - 10. * nb_s_e1);
    size_t nb_mill = int((t_time - int(t_time)) * 1000. + 0.499999);
    size_t nb_mill_e2 = int(nb_mill / 100.);
    size_t nb_mill_e1 = int((nb_mill - 100 * nb_mill_e2) / 10.);
    size_t nb_mill_e0 = int(nb_mill - 100 * nb_mill_e2 - 10 * nb_mill_e1);

    file_name += std::to_string(nb_d_e3)
              + std::to_string(nb_d_e2) + std::to_string(nb_d_e1)
              + std::to_string(nb_d_e0) + "d"
              + std::to_string(nb_h_e1) + std::to_string(nb_h_e0) + "h"
              + std::to_string(nb_m_e1) + std::to_string(nb_m_e0) + "m"
              + std::to_string(nb_s_e1) + std::to_string(nb_s_e0) + "s"
              + std::to_string(nb_mill_e2) + std::to_string(nb_mill_e1)
              + std::to_string(nb_mill_e0) + label + ".txt";
    file_name = absolutePath(file_name);
}


void swIO::createEmptyFolders(const std::string &log_file)
/**
 * @brief creates log and output folder
 * @param[in] log_file: (meant to be part of log files), could be any log file (converter, solver, etc.)
 */
{
    createEmptyLogFolder(log_file);
    createEmptyOutputFolder(false);
}


void swIO::createEmptyLogFolder(const std::string &log_file)
/**
 * @brief Creates the log folder (and removes possible existing file)
 * @param[in] log_file (meant to be part of log files), could be any log file (converter, solver, etc.)
 */
{
    std::string log_path = absolutePath(swFile::log_folder);

    if (std::filesystem::exists(log_path) == true)   // the 'log' folder already exists
    {
        std::string log_file_path = absolutePath(log_file);

        if (std::filesystem::exists(log_file_path) == true)   // the 'file_type' file already exists => remove
        {
            std::filesystem::remove(log_file_path);
        }
    }
    else
    {
        // the 'log' folder does not exists => create
        std::filesystem::create_directory(log_path);
    }
}


void swIO::createEmptyOutputFolder(bool remove_output)
/**
 * @brief creates output folder 
 * @param[in] remove_output to decide if existing files should be removed from output folder
 * 
 */
{
    //---- output folder
    std::string output_path = absolutePath(swFile::output_folder);
    if (remove_output == true)
    {
        std::filesystem::remove_all(output_path);
    }
    std::filesystem::create_directory(output_path);
    //==== output folder
}


bool swIO::createTxtFileMapDecimal(std::string const& file_name, std::vector<std::string> const& var_names,
    std::vector<std::vector<double>> const& var_values, size_t precision)
/**
 * @brief Creates and writes a map file containing float/double vector values
 * @param[in] file_name: name of the txt file in which results should be stored
 * @param[in] var_names: name of the outputs to be stored in the header
 * @param[in] var_values: numerical values for the outputs to be stored after the header (one line per cell)
 * @return EXIT_SUCCESS / FAILURE
*/
{
    std::ofstream file_w;
    file_w.open(file_name);
    //-- Header
    size_t nb_char = precision + 7; 
    for (size_t i = 0; i < var_names.size(); i++)
    {
        file_w << var_names[i];
        if (var_names[i].size() < nb_char)
        {
            for (size_t k = 0; k < nb_char - var_names[i].size(); k++)
            {
                file_w << " ";
            }
        }
        file_w << "\t";
    }
    file_w << "\n";
    //== Header
    
    //-- Values
    for (size_t idx_record = 0; idx_record < var_values[0].size(); idx_record++)
    {
        for (size_t idx_variable = 0; idx_variable < var_values.size(); idx_variable++)
        {
            file_w << swString::numberToStr(var_values[idx_variable][idx_record], 'e', precision)  << "\t";
        }
        file_w << "\n";
    }
    file_w.close();
    //== Values
    return EXIT_SUCCESS;
}


bool swIO::createFileMapString(std::string const& file_name, std::vector<std::string> const& var_names,
    std::vector<std::vector<std::string>> const& var_values)
/**
 * @brief Creates and writes a map file containing string values
 * @param[in] file_name The name of the file to be created
 * @param[in] var_names[i]: the name of the i-th variable
 * @param[in] var_values[i][n] the value of the n-th record for the i-th variable, cast to string to accommodate for any type
 * @return EXIT_SUCCESS / FAILURE
*/
{
    std::ofstream file_w;
    file_w.open(file_name);
    //-- Header
    for (size_t i = 0; i < var_names.size(); i++)
    {
        file_w << var_names[i] << "\t";
    }
    file_w << "\n";
    //== Header
     //-- Values
    for (size_t idx_record = 0; idx_record < var_values[0].size(); idx_record++)
    {
        for (size_t idx_variable = 0; idx_variable < var_values.size(); idx_variable++)
        {
            file_w <<  var_values[idx_variable][idx_record] << "\t";
        }
        file_w << "\n";
    }
    file_w.close();
    //== Values
    return EXIT_SUCCESS;
}


void swIO::errorReadingLine(std::ifstream& file,
                            const std::string& filename,
                            const std::string& line_txt,
                            int line_input_file, int line_code,
                            const std::string& function_name,
                            const std::string& cpp_name)
/**
 * @brief Output detailed message if error while reading a line of the input file + fatal error
 * @param[in] file: stream handler of the file
 * @param[in] filename: the name of the file that is being read
 * @param[in] line_txt: the line that could not be correctly read
 * @param[in] line_input_file: line number of the input file in which the error ocurred
 * @param[in] line_code: line of the reading function in which the error occured (use __LINE__)
 * @param[in] function_name: name of the function in which the error occurred (use __func__)
 * @param[in] cpp_name: name of the cpp in which the calling function is declared (use __FILE__)
 */
{
    file.close();
    logError() << "Can't read file " << filename.c_str()
               << ", line=" << line_input_file
               << ": " << line_txt.c_str();
    logFatal(cpp_name, function_name.c_str(), line_code);
}


bool swIO::fillBinFileMapDecimal(const std::string &filename, double time, 
                            std::vector<std::vector<double>> const& var_values)
/**
 * @brief store results in the dedicated binary file (one file for all time steps)
 * @param[in] filename: name of the output binary file
 * @param[in] time: time at which information is stored
 * @param[in] var_values: table for all variable values, for every cell
 * @return EXIT_SUCCESS / FAILURE
 */
{
    std::ofstream binary_file_stream;
    binary_file_stream.open(filename, std::ios::app | std::ios::binary);
    binary_file_stream.write(reinterpret_cast<const char*>(&time), sizeof(time));
    for (size_t j = 0; j < var_values[0].size(); j++)
    {
        for (size_t i = 2 ; i < var_values.size() ; i++) // Do not write x and y (start at i=2)
        {
            binary_file_stream.write(reinterpret_cast<const char*>(&var_values[i][j]), sizeof(var_values[i][j]));
        }
    }
    binary_file_stream.close();

    return EXIT_SUCCESS;
}


bool swIO::initBinFileMapDecimal(const std::string &filename, const std::string &geo_filename,
                                 const std::vector<std::string> &var_names,
                                 const double t_0, const double t_max, const double dt_map)
/**
 * @brief opens the binary result file and store preliminary info (geometry content, variable names, initial and final times...)
 * @param[in] filename: name of the output binary file
 * @param[in] geo_filename: name of the geometry file
 * @param[in] var_names: table for all variable names
 * @param[in] t_0: initial time
 * @param[in] t_max: final time
 * @param[in] dt_map: time interval between two consecutive storages
 * @return EXIT_SUCCESS / FAILURE
*/ 
{
        //--- Create the file and open the stream
        std::ofstream binary_file_stream;
        binary_file_stream.open(filename, std::ios::binary);

        if (binary_file_stream.is_open() == false)
        {
            logError() << filename << "cannot be opened";
            logFatal(__FILE__, __func__, __LINE__);
        }
        //=== Create the file and open the stream

        //--- Write geometry
        // Open geometry.sw2d
        std::ifstream geo_stream(geo_filename, std::ios::binary);
        if (geo_stream.is_open() == false)
        {
            logError() << "can't open file " << geo_filename;
            logFatal(__FILE__, __func__, __LINE__);
        }

        // copy version number
        size_t version_lenght;
        geo_stream.read(reinterpret_cast<char*>(&version_lenght), sizeof(version_lenght));
        std::string version(version_lenght, '\0');
        geo_stream.read(&version[0], version_lenght);
        binary_file_stream.write(reinterpret_cast<const char*>(&version_lenght), sizeof(version_lenght));
        binary_file_stream.write(version.c_str(), version_lenght);
        
        //-- read number of elements
        size_t nb_cells;
        size_t nb_edges_boundary;
        size_t nb_edges_intern;
        size_t nb_nodes;
        size_t nb_boundaries;

        geo_stream.read(reinterpret_cast<char*>(&nb_cells), sizeof(nb_cells));
        geo_stream.read(reinterpret_cast<char*>(&nb_edges_boundary), sizeof(nb_edges_boundary));
        geo_stream.read(reinterpret_cast<char*>(&nb_edges_intern), sizeof(nb_edges_intern));
        geo_stream.read(reinterpret_cast<char*>(&nb_nodes), sizeof(nb_nodes));
        geo_stream.read(reinterpret_cast<char*>(&nb_boundaries), sizeof(nb_boundaries));
        binary_file_stream.write(reinterpret_cast<const char*>(&nb_cells), sizeof(nb_cells));
        binary_file_stream.write(reinterpret_cast<const char*>(&nb_edges_boundary), sizeof(nb_edges_boundary));
        binary_file_stream.write(reinterpret_cast<const char*>(&nb_edges_intern), sizeof(nb_edges_intern));
        binary_file_stream.write(reinterpret_cast<const char*>(&nb_nodes), sizeof(nb_nodes));
        binary_file_stream.write(reinterpret_cast<const char*>(&nb_boundaries), sizeof(nb_boundaries));
        //== read number of elements
        
        //-- read cells
        std::vector<size_t> nb_edges_bnd_per_cell;
        std::vector<size_t> nb_edges_int_per_cell;
        std::vector<size_t> nb_nodes_per_cell;
        std::vector<size_t> nb_neigh_per_cell;

        for (size_t n = 0; n < nb_cells; ++n)
        {      
            size_t idx;
            size_t sms_id;
            size_t nb_c_node;
            size_t nb_c_edge_boundary;
            size_t nb_c_edge_intern;
            size_t nb_neigh;
            double area;
            double perimeter;
            double xg;
            double yg;
            double zb;

            geo_stream.read(reinterpret_cast<char*>(&idx), sizeof(idx));
            geo_stream.read(reinterpret_cast<char*>(&sms_id), sizeof(sms_id));
            geo_stream.read(reinterpret_cast<char*>(&nb_c_node), sizeof(nb_c_node));
            geo_stream.read(reinterpret_cast<char*>(&nb_c_edge_boundary), sizeof(nb_c_edge_boundary));
            geo_stream.read(reinterpret_cast<char*>(&nb_c_edge_intern), sizeof(nb_c_edge_intern));
            geo_stream.read(reinterpret_cast<char*>(&nb_neigh), sizeof(nb_neigh));
            geo_stream.read(reinterpret_cast<char*>(&area), sizeof(area));
            geo_stream.read(reinterpret_cast<char*>(&perimeter), sizeof(perimeter));
            geo_stream.read(reinterpret_cast<char*>(&xg), sizeof(xg));
            geo_stream.read(reinterpret_cast<char*>(&yg), sizeof(yg));
            geo_stream.read(reinterpret_cast<char*>(&zb), sizeof(zb));

            binary_file_stream.write(reinterpret_cast<const char*>(&idx), sizeof(idx));
            binary_file_stream.write(reinterpret_cast<const char*>(&sms_id), sizeof(sms_id));
            binary_file_stream.write(reinterpret_cast<const char*>(&nb_c_node), sizeof(nb_c_node));
            binary_file_stream.write(reinterpret_cast<const char*>(&nb_c_edge_boundary), sizeof(nb_c_edge_boundary));
            binary_file_stream.write(reinterpret_cast<const char*>(&nb_c_edge_intern), sizeof(nb_c_edge_intern));
            binary_file_stream.write(reinterpret_cast<const char*>(&nb_neigh), sizeof(nb_neigh));
            binary_file_stream.write(reinterpret_cast<const char*>(&area), sizeof(area));
            binary_file_stream.write(reinterpret_cast<const char*>(&perimeter), sizeof(perimeter));
            binary_file_stream.write(reinterpret_cast<const char*>(&xg), sizeof(xg));
            binary_file_stream.write(reinterpret_cast<const char*>(&yg), sizeof(yg));
            binary_file_stream.write(reinterpret_cast<const char*>(&zb), sizeof(zb));

            nb_edges_bnd_per_cell.push_back(nb_c_edge_boundary);
            nb_edges_int_per_cell.push_back(nb_c_edge_intern);
            nb_nodes_per_cell.push_back(nb_c_node);
            nb_neigh_per_cell.push_back(nb_neigh);
        }
        //== read cells

        //-- read nodes
        std::vector<size_t> nb_edges_per_point;
        std::vector<size_t> nb_cells_per_point;

        for (size_t p = 0; p < nb_nodes; ++p)
        {        
            size_t idx;
            size_t sms_id;
            size_t ncelp;
            size_t nedgp;
            double xnod;
            double ynod;
            double znod;

            geo_stream.read(reinterpret_cast<char*>(&idx), sizeof(idx));
            geo_stream.read(reinterpret_cast<char*>(&sms_id), sizeof(sms_id));
            geo_stream.read(reinterpret_cast<char*>(&xnod), sizeof(xnod));
            geo_stream.read(reinterpret_cast<char*>(&ynod), sizeof(ynod));
            geo_stream.read(reinterpret_cast<char*>(&znod), sizeof(znod));
            geo_stream.read(reinterpret_cast<char*>(&ncelp), sizeof(ncelp));
            geo_stream.read(reinterpret_cast<char*>(&nedgp), sizeof(nedgp));

            binary_file_stream.write(reinterpret_cast<const char*>(&idx), sizeof(idx));
            binary_file_stream.write(reinterpret_cast<const char*>(&sms_id), sizeof(sms_id));
            binary_file_stream.write(reinterpret_cast<const char*>(&xnod), sizeof(xnod));
            binary_file_stream.write(reinterpret_cast<const char*>(&ynod), sizeof(ynod));
            binary_file_stream.write(reinterpret_cast<const char*>(&znod), sizeof(znod));
            binary_file_stream.write(reinterpret_cast<const char*>(&ncelp), sizeof(ncelp));
            binary_file_stream.write(reinterpret_cast<const char*>(&nedgp), sizeof(nedgp));

            nb_edges_per_point.push_back(nedgp);
            nb_cells_per_point.push_back(ncelp);
        }
        //== read nodes

        //-- read boundary edges
        for (size_t i = 0; i < nb_edges_boundary; ++i)
        {
            size_t idx;
            size_t node_beg;
            size_t node_end;
            double x_int;
            double y_int;
            double z_int;
            double d_norm;
            double r_lenght;
            double norm_x;
            double norm_y;
            int bc_code;
        
            geo_stream.read(reinterpret_cast<char*>(&idx), sizeof(idx));
            geo_stream.read(reinterpret_cast<char*>(&node_beg), sizeof(node_beg));
            geo_stream.read(reinterpret_cast<char*>(&node_end), sizeof(node_end));
            geo_stream.read(reinterpret_cast<char*>(&x_int), sizeof(x_int));
            geo_stream.read(reinterpret_cast<char*>(&y_int), sizeof(y_int));
            geo_stream.read(reinterpret_cast<char*>(&z_int), sizeof(z_int));
            geo_stream.read(reinterpret_cast<char*>(&d_norm), sizeof(d_norm));
            geo_stream.read(reinterpret_cast<char*>(&r_lenght), sizeof(r_lenght));
            geo_stream.read(reinterpret_cast<char*>(&norm_x), sizeof(norm_x));
            geo_stream.read(reinterpret_cast<char*>(&norm_y), sizeof(norm_y));
            geo_stream.read(reinterpret_cast<char*>(&bc_code), sizeof(bc_code));

            binary_file_stream.write(reinterpret_cast<char*>(&idx), sizeof(idx));
            binary_file_stream.write(reinterpret_cast<char*>(&node_beg), sizeof(node_beg));
            binary_file_stream.write(reinterpret_cast<char*>(&node_end), sizeof(node_end));
            binary_file_stream.write(reinterpret_cast<char*>(&x_int), sizeof(x_int));
            binary_file_stream.write(reinterpret_cast<char*>(&y_int), sizeof(y_int));
            binary_file_stream.write(reinterpret_cast<char*>(&z_int), sizeof(z_int));
            binary_file_stream.write(reinterpret_cast<char*>(&d_norm), sizeof(d_norm));
            binary_file_stream.write(reinterpret_cast<char*>(&r_lenght), sizeof(r_lenght));
            binary_file_stream.write(reinterpret_cast<char*>(&norm_x), sizeof(norm_x));
            binary_file_stream.write(reinterpret_cast<char*>(&norm_y), sizeof(norm_y));
            binary_file_stream.write(reinterpret_cast<char*>(&bc_code), sizeof(bc_code));

            size_t cell_idx;
            geo_stream.read(reinterpret_cast<char*>(&cell_idx), sizeof(cell_idx));
            binary_file_stream.write(reinterpret_cast<char*>(&cell_idx), sizeof(cell_idx)); 
        }
        //== read boundary edges

        //-- read internal edges
        for (size_t i = 0; i < nb_edges_intern; ++i)
        {
            size_t idx;
            size_t node_beg;
            size_t node_end;
            double x_int;
            double y_int;
            double z_int;
            double r_lenght;
            double d_norm;
            double norm_x;
            double norm_y;
        
            geo_stream.read(reinterpret_cast<char*>(&idx), sizeof(idx));
            geo_stream.read(reinterpret_cast<char*>(&node_beg), sizeof(node_beg));
            geo_stream.read(reinterpret_cast<char*>(&node_end), sizeof(node_end));
            geo_stream.read(reinterpret_cast<char*>(&x_int), sizeof(x_int));
            geo_stream.read(reinterpret_cast<char*>(&y_int), sizeof(y_int));
            geo_stream.read(reinterpret_cast<char*>(&z_int), sizeof(z_int));
            geo_stream.read(reinterpret_cast<char*>(&d_norm), sizeof(d_norm));
            geo_stream.read(reinterpret_cast<char*>(&r_lenght), sizeof(r_lenght));
            geo_stream.read(reinterpret_cast<char*>(&norm_x), sizeof(norm_x));
            geo_stream.read(reinterpret_cast<char*>(&norm_y), sizeof(norm_y));

            binary_file_stream.write(reinterpret_cast<char*>(&idx), sizeof(idx));
            binary_file_stream.write(reinterpret_cast<char*>(&node_beg), sizeof(node_beg));
            binary_file_stream.write(reinterpret_cast<char*>(&node_end), sizeof(node_end));
            binary_file_stream.write(reinterpret_cast<char*>(&x_int), sizeof(x_int));
            binary_file_stream.write(reinterpret_cast<char*>(&y_int), sizeof(y_int));
            binary_file_stream.write(reinterpret_cast<char*>(&z_int), sizeof(z_int));
            binary_file_stream.write(reinterpret_cast<char*>(&d_norm), sizeof(d_norm));
            binary_file_stream.write(reinterpret_cast<char*>(&r_lenght), sizeof(r_lenght));
            binary_file_stream.write(reinterpret_cast<char*>(&norm_x), sizeof(norm_x));
            binary_file_stream.write(reinterpret_cast<char*>(&norm_y), sizeof(norm_y));

            for (size_t n = 0; n < 2 ; ++n)
            {
                size_t cell_idx;
                geo_stream.read(reinterpret_cast<char*>(&cell_idx), sizeof(cell_idx));
                binary_file_stream.write(reinterpret_cast<char*>(&cell_idx), sizeof(cell_idx));
            }
        }
        //== read internal edges

        //-- read boundaries
        for (size_t k = 0; k < nb_boundaries; k++)
        {
            size_t bnd_size;
            geo_stream.read(reinterpret_cast<char*>(&bnd_size), sizeof(bnd_size));
            binary_file_stream.write(reinterpret_cast<char*>(&bnd_size), sizeof(bnd_size));

            for (size_t i = 0; i < bnd_size; i++)
            {
                size_t edge_idx;
                geo_stream.read(reinterpret_cast<char*>(&edge_idx), sizeof(edge_idx));
                binary_file_stream.write(reinterpret_cast<char*>(&edge_idx), sizeof(edge_idx));
            }
        }

        for (size_t k = 0; k < nb_boundaries; k++)
        {
            int code;
            geo_stream.read(reinterpret_cast<char*>(&code), sizeof(code));
            binary_file_stream.write(reinterpret_cast<char*>(&code), sizeof(code));
        }
        //== read boundaries

        //-- read cell link with other elements
        for (size_t n = 0; n < nb_cells; n++)
        {
            for (size_t p = 0; p < nb_nodes_per_cell[n]; ++p)
            {
                size_t node_pt;
                geo_stream.read(reinterpret_cast<char*>(&node_pt), sizeof(node_pt));
                binary_file_stream.write(reinterpret_cast<const char*>(&node_pt), sizeof(node_pt));
            }

            for (size_t i = 0; i < nb_edges_bnd_per_cell[n]; ++i)
            {
                size_t edge_bnd_pt;
                geo_stream.read(reinterpret_cast<char*>(&edge_bnd_pt), sizeof(edge_bnd_pt));
                binary_file_stream.write(reinterpret_cast<const char*>(&edge_bnd_pt), sizeof(edge_bnd_pt));
            }

            for (size_t i = 0; i < nb_edges_int_per_cell[n]; ++i)
            {
                size_t edge_int_pt;
                geo_stream.read(reinterpret_cast<char*>(&edge_int_pt), sizeof(edge_int_pt));
                binary_file_stream.write(reinterpret_cast<const char*>(&edge_int_pt), sizeof(edge_int_pt));
            }

            for (size_t n_n = 0; n_n < nb_neigh_per_cell[n]; ++n_n)
            {
                size_t cell_pt;
                geo_stream.read(reinterpret_cast<char*>(&cell_pt), sizeof(cell_pt));
                binary_file_stream.write(reinterpret_cast<const char*>(&cell_pt), sizeof(cell_pt));
            }
        }
        //== read cell link with other elements

        //-- read node link with other elements
        for (size_t p=0; p<nb_nodes; p++)
        {
            for (size_t n_p = 0; n_p < nb_cells_per_point[p]; ++n_p)
            {
                size_t cel_pt;
                geo_stream.read(reinterpret_cast<char*>(&cel_pt), sizeof(cel_pt));
                binary_file_stream.write(reinterpret_cast<const char*>(&cel_pt), sizeof(cel_pt));
            }
        }
        //== read node link with other elements
        //--- Read the geometry
        //=== Write geometry
        
        //--- Write variable names
        size_t nb_var = var_names.size() + 1;
        binary_file_stream.write(reinterpret_cast<const char*>(&nb_var), sizeof(nb_var));

        std::string zb_name = "bottom_elevation";
        size_t zb_name_lenght = zb_name.size();
        binary_file_stream.write(reinterpret_cast<const char*>(&zb_name_lenght), sizeof(zb_name_lenght));
        binary_file_stream.write(zb_name.c_str(), zb_name_lenght);
        
        for (std::string name : var_names)
        {
            size_t name_lenght = name.size();
            binary_file_stream.write(reinterpret_cast<const char*>(&name_lenght), sizeof(name_lenght));
            binary_file_stream.write(name.c_str(), name_lenght);
        }
        //=== Write variable names

        //--- Write time information
        binary_file_stream.write(reinterpret_cast<const char*>(&t_0), sizeof(t_0));
        binary_file_stream.write(reinterpret_cast<const char*>(&t_max), sizeof(t_max));
        binary_file_stream.write(reinterpret_cast<const char*>(&dt_map), sizeof(dt_map));
        //=== Write time information

        binary_file_stream.close();
        
    return EXIT_SUCCESS;
}


bool swIO::readBoundaryRatingCurveFile(const std::string& filename, std::vector<std::vector<std::string>>& rating_curve_data)
/**
 * @brief reads rating curve data from dedicated input file
 * @param[in] filename: path to the hydro_boundary_rating_curve file
 * @param[in,out] rating_curve_data: string content of the input file
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    std::ifstream file_bc;
    file_bc.open(filename);

    if (file_bc.is_open() == false)
    {
        logError() << filename << "cannot be opened";
        logFatal(__FILE__, __func__, __LINE__);
    }

    size_t bc_nb = 0;
    std::string line_txt;

    /* Read 1 commented lines */
    line_txt = swIO::readLine(file_bc);

    while (file_bc.eof() == false)
    {
        rating_curve_data.resize(bc_nb+1);
        /* Read 1 commented lines */
        line_txt = swIO::readLine(file_bc);
        // Read 5 items: Code, type, z0, line1, line2
        for (size_t i = 0; i < 5 ; i++)
        {
            line_txt = swIO::readLine(file_bc);
            rating_curve_data[bc_nb].push_back(line_txt);
        }
        bc_nb++;
    }
    return EXIT_SUCCESS;
}


bool swIO::readBoundaryTimeSeriesFile(const std::string& filename,
                                      std::vector<int> &bc_timeseries_codes,
                                      std::vector<std::string> &bc_timeseries_type,
                                      std::vector<double> &timebc,
                                      std::vector<std::vector<double>> &val_bc)
/**
 * @brief swIO::readBoundaryTimeSeriesFile
 * @param[in] filename: path to the hydro_boundary_timeseries file
 * @param[in] bc_timeseries_list: list of the timeseries bc
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    std::ifstream file_bc;
    file_bc.open(filename);

    if (file_bc.is_open() == false)
    {
        logError() << filename << "cannot be opened";
        logFatal(__FILE__, __func__, __LINE__);
    }

    size_t line_nb = 0;
    std::string line_txt;
    bool conversion_success = true;

    /* Read 2 commented lines */

    line_txt = swIO::readLine(file_bc);
    ++line_nb;
    line_txt = swIO::readLine(file_bc);
    ++line_nb;

    /* Read the number of BC types */
    line_txt = swIO::readLine(file_bc);
    ++line_nb;
    size_t nseries = swString::stringToUInt(swString::trim(line_txt), conversion_success);
    if (conversion_success == false)
    {
        logError() << "line is misformatted: format conversion";
        swIO::errorReadingLine(file_bc, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
    }

    /* Read 1 commented line */
    line_txt = swIO::readLine(file_bc);
    ++line_nb;

    /* Read the bc_codes */
    line_txt = swIO::readLine(file_bc);
    ++line_nb;
    const std::vector<std::string>& bc_codes_string = swString::splitStr(swString::trim(line_txt), swString::regexWhiteSpace());
    if (bc_codes_string.size() != nseries)
    {
        logError() << "incorrect number of specified boundary conditions";
        logError() << nseries << "defined at line 3";
        logError() << bc_codes_string.size() << "bc_codes defined at line 5";
        logFatal(__FILE__, __func__, __LINE__);
    }
    for (size_t k = 0; k < nseries; ++k)
    {
        int code = swString::stringToInt(bc_codes_string[k], conversion_success);
        if (conversion_success == false)
        {
            logError() << "line is misformatted: format conversion";
            swIO::errorReadingLine(file_bc, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
        }
        bc_timeseries_codes.push_back(code);
    }

    /* Read 1 commented line */
    line_txt = swIO::readLine(file_bc);
    ++line_nb;

    /* Read the bc type */
    line_txt = swIO::readLine(file_bc);
    ++line_nb;
    const std::vector<std::string>& bc_types_string = swString::splitStr(swString::trim(line_txt), swString::regexWhiteSpace());
    if (bc_types_string.size() != nseries)
    {
        logError() << "incorrect number of specified boundary conditions";
        logError() << nseries << "defined at line 3";
        logError() << bc_types_string.size() << "bc_types defined at line 7";
        logFatal(__FILE__, __func__, __LINE__);
    }
    for (size_t k = 0; k < nseries; ++k)
    {
        bc_timeseries_type.push_back(bc_types_string[k]);
    }

    /* Read 2 commented lines */
    line_txt = swIO::readLine(file_bc);
    ++line_nb;
    line_txt = swIO::readLine(file_bc);
    ++line_nb;

    if (file_bc.eof() == true)
    {
        logError() << "The file " << filename << " does not contain any boundary condition values";
        logFatal(__FILE__, __func__, __LINE__);
    }

    //// For every type, read the data
    val_bc.resize(nseries);
    while (file_bc.eof() == false)
    {
        line_txt = swIO::readLine(file_bc);
        ++line_nb;

        if (line_txt.size() != 0) // Avoid empty lines
        {
            std::vector<std::string> vals_str = swString::splitStr(swString::trim(line_txt), swString::regexWhiteSpace());
            double t = swString::stringToDouble(vals_str[0], conversion_success);
            if (conversion_success == false)
            {
                logError() << "line is misformatted: format conversion";
                swIO::errorReadingLine(file_bc, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
            }
            timebc.push_back(t);
            for (size_t j = 0; j < nseries; ++j)
            {
                double val = swString::stringToDouble(vals_str[j + 1], conversion_success);
                if (conversion_success == false)
                {
                    logError() << "line is misformatted: format conversion";
                    swIO::errorReadingLine(file_bc, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
                }
                val_bc[j].push_back(val);
            }
        }
    }
    return EXIT_SUCCESS;
}


bool swIO::readTimeSeriesScalarFile(const std::string& filename,
                                    std::vector<double>& time_vector,
                                    std::vector<std::vector<double>>& val_table)
/**
 * @brief Reads a couple of tabulated (time_vector, scalar val_table)
 * @param[in] The name of the file that contains the time sries
 * @param[in] time_vector[i]: the tabulated time for record number i
 * @param[in] val_table [n][i]: the value for the n-th scalar time series and record i
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    std::ifstream file_series;
    file_series.open(filename);

    if (file_series.is_open() == false)
    {
        logError() << filename << "cannot be opened";
        logFatal(__FILE__, __func__, __LINE__);
    }

    size_t line_nb = 0;
    std::string line_txt;
    bool conversion_success = true;

    /* Read 2 commented lines */

    line_txt = swIO::readLine(file_series);
    ++line_nb;
    line_txt = swIO::readLine(file_series);
    ++line_nb;

    /* Read the number of time series */
    line_txt = swIO::readLine(file_series);
    ++line_nb;
    size_t nseries = swString::stringToUInt(swString::trim(line_txt), conversion_success);
    if (conversion_success == false)
    {
        logError() << "line is misformatted: format conversion";
        swIO::errorReadingLine(file_series, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
    }

    /* Read 2 commented lines */
    line_txt = swIO::readLine(file_series);
    ++line_nb;
    line_txt = swIO::readLine(file_series);
    ++line_nb;

    if (file_series.eof() == true)
    {
        logError() << "The file " << filename << " does not contain any boundary condition values";
        logFatal(__FILE__, __func__, __LINE__);
    }

    //// Read data
    val_table.resize(nseries);
    while (file_series.eof() == false)
    {
        line_txt = swIO::readLine(file_series);
        ++line_nb;

        if (line_txt.size() != 0) // Avoid empty lines
        {
            std::vector<std::string> vals_str = swString::splitStr(swString::trim(line_txt), swString::regexWhiteSpace());
            double t = swString::stringToDouble(vals_str[0], conversion_success);
            if (conversion_success == false)
            {
                logError() << "line is misformatted: format conversion";
                swIO::errorReadingLine(file_series, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
            }
            time_vector.push_back(t);
            for (size_t j = 0; j < nseries; ++j)
            {
                double val = swString::stringToDouble(vals_str[j + 1], conversion_success);
                if (conversion_success == false)
                {
                    logError() << "line is misformatted: format conversion";
                    swIO::errorReadingLine(file_series, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
                }
                val_table[j].push_back(val);
            }
        }
    }
    return EXIT_SUCCESS;
}


bool swIO::readTimeSeries2dVectorFile(const std::string& filename,
                                      std::vector<double>& time_vector,
                                      std::vector<std::vector<double>>& val_1_table,
                                      std::vector<std::vector<double>>& val_2_table)
/**
 * @brief Reads a triplet of tabulated (time_vector, val_1_table, val_2_table) - val_1 and val_2 may be the 2 components of a vector
 * @param[in] The name of the file that contains the time sries
 * @param[out] time_vector[n]: the tabulated time for record number n
 * @param[out] val_1_table [i][n]: the value of the i-th val_1 variable for the n-th record
 * @param[out] val_2_table [i][n]: the value of the i-th val_1 variable for the n-th record
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    std::ifstream file_series;
    file_series.open(filename);

    if (file_series.is_open() == false)
    {
        logError() << filename << "cannot be opened";
        logFatal(__FILE__, __func__, __LINE__);
    }

    size_t line_nb = 0;
    std::string line_txt;
    bool conversion_success = true;

    /* Read 2 commented lines */

    line_txt = swIO::readLine(file_series);
    ++line_nb;
    line_txt = swIO::readLine(file_series);
    ++line_nb;

    /* Read the number of time series */
    line_txt = swIO::readLine(file_series);
    ++line_nb;
    size_t nseries = swString::stringToUInt(swString::trim(line_txt), conversion_success);
    if (conversion_success == false)
    {
        logError() << "line is misformatted: format conversion";
        swIO::errorReadingLine(file_series, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
    }

    /* Read 2 commented lines */
    line_txt = swIO::readLine(file_series);
    ++line_nb;
    line_txt = swIO::readLine(file_series);
    ++line_nb;

    if (file_series.eof() == true)
    {
        logError() << "The file " << filename << " does not contain any boundary condition values";
        logFatal(__FILE__, __func__, __LINE__);
    }

    //// Read data
    val_1_table.resize(nseries);
    val_2_table.resize(nseries);
    while (file_series.eof() == false)
    {
        line_txt = swIO::readLine(file_series);
        ++line_nb;

        if (line_txt.size() != 0) // Avoid empty lines
        {
            std::vector<std::string> vals_str = swString::splitStr(swString::trim(line_txt), swString::regexWhiteSpace());
            double t = swString::stringToDouble(vals_str[0], conversion_success);
            if (conversion_success == false)
            {
                logError() << "line is misformatted: format conversion";
                swIO::errorReadingLine(file_series, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
            }
            time_vector.push_back(t);
            for (size_t j = 0; j < nseries; ++j)
            {
                // Convert 1rst component of vector
                double val = swString::stringToDouble(vals_str[2 * j + 1], conversion_success);
                if (conversion_success == false)
                {
                    logError() << "line is misformatted: format conversion";
                    swIO::errorReadingLine(file_series, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
                }
                val_1_table[j].push_back(val);
                // Convert 2nd component of vector
                val = swString::stringToDouble(vals_str[2 * j + 2], conversion_success);
                if (conversion_success == false)
                {
                    logError() << "line is misformatted: format conversion";
                    swIO::errorReadingLine(file_series, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
                }

                val_2_table[j].push_back(val);
            }
        }
    }
    return EXIT_SUCCESS;
}


void swIO::readVectorParameterFile(std::string const& file_name,
                                          size_t const& n_params, size_t const& n_records,
                                          std::vector<std::string>& param_name,
                                          std::vector<std::vector<double>>& param_array)
/**
 * @brief Reads an input file containing a known number of parameters for a known number of records
 * @param[in] file_name : name of the file to be read
 * @param[in] n_params : the number of parameters (columns) in the file
 * @param[in] n_records : the number of records to be filled (number of lines in the file if parameter is not uniform)
 * @param[in] param_array [n][k] : the numerical value of the k-th record of the n-th parameter in the file
 */
{
    //---- read file
    std::ifstream file_r;
    file_r.open(file_name);
    if (file_r.is_open() == false)
    {
        logError() << "    Unable to open file " << file_name;
        logFatal(__FILE__, __func__, __LINE__);
    }

    logInfo() << "Importing parameters from file:" << file_name;

    size_t line_nb = 0;
    std::string line_txt;

    //-- Prepare arrays
    param_array.resize(n_params);
    param_name.resize(n_params);
    //== Prepare arrays

    //------ Read distribution type (uniform or not)
    line_txt = swIO::readLine(file_r);
    ++line_nb;
    std::vector<std::string> line_split = swIO::stringSplit(line_txt);
    if (line_split.size() != 2)
    {
        logError() << "line " << line_nb << " is misformatted: number of arguments = " << line_split.size();
        logFatal(__FILE__, __func__, __LINE__);
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
        logError() << "Incorrect value for 'Unif' parameter";
        logError() << "  allowed value for parameter 'Unif' : 0 or 1";
        logFatal(__FILE__, __func__, __LINE__);
    }
    logInfo() << "Parameter uniformity = " << unif;
    //====== Read distribution type (uniform or not)


    //------ Read dummy line ("==== Default Param"
    line_txt = swIO::readLine(file_r);
    ++line_nb;
    //====== Read dummy line ("==== Default Param"

    //------ Read header (name of variables)
    line_txt = swIO::readLine(file_r);
    ++line_nb;
    line_split = swIO::stringSplit(line_txt);
    if (line_split.size() != n_params)
    {
        logError() << "line " << line_nb << " is misformatted";
        logError() << "  expecting" << n_params << "parameters (" << line_split.size() << "found)";
        logFatal(__FILE__, __func__, __LINE__);
    }
    //-- initialize vector for storing header
    for (size_t i = 0; i < line_split.size(); i++)
    {
        param_name[i] = line_split[i];
        logInfo() << "  param number and name = " << i << ", " << param_name[i];
    }
    //== initialize vector for storing header
    //======= Read header (name of variables)


    std::vector<double> vals_line = std::vector<double>(n_params);

    //------ Uniform case
    if (unif == 1)
    {
        line_txt = swIO::readLine(file_r);
        ++line_nb;
        line_split = swIO::stringSplit(line_txt);
        if (line_split.size() != n_params)
        {
            logError() << "line " << line_nb << " is misformatted";
            logError() << "  expecting" << n_params << "parameters (" << line_split.size() << "found)";
            logFatal(__FILE__, __func__, __LINE__);
        }
        else
        {
            for (size_t idx_param = 0; idx_param < n_params; ++idx_param)
            {
                param_array[idx_param] = { };
                bool conversion_success = true;
                double val = swString::stringToDouble(line_split[idx_param], conversion_success);
                logInfo() << "Assigning uniform variable - record, value = " << idx_param << val;
                if (conversion_success == false)
                {
                    logError() << "line is misformatted: format conversion";
                    swIO::errorReadingLine(file_r, file_name, line_txt, line_nb, __LINE__, __func__, __FILE__);
                }
                for (size_t idx_record = 0; idx_record < n_records; idx_record++)
                {
                    param_array[idx_param].push_back(val);
                }
            }
        }
        file_r.close();
    }
    //====== Uniform case
    // 
    // 
    //------ Variable case
    else
    {
        /* Read 2 commented lines */
        line_txt = swIO::readLine(file_r);// Default values
        ++line_nb;
        line_txt = swIO::readLine(file_r);// ==== Distrib
        ++line_nb;

        const size_t nb_dummy_lines = 5;
        for (size_t idx_param = 0; idx_param < n_params; idx_param++)
        {
            param_array[idx_param] = { };
        }
        for (size_t idx_record = 0; idx_record < n_records; idx_record++)
        {
            if (file_r.eof() == true)
            {
                file_r.close();
                logError() << "Not enough values in file" << file_name;
                logError() << "Expected/existing number of lines containing variable parameters:"
                    << n_records << line_nb - nb_dummy_lines - 1;
                logFatal(__FILE__, __func__, __LINE__);
            }

            line_txt = swIO::readLine(file_r);
            ++line_nb;
            line_split = swIO::stringSplit(line_txt);
            if (line_split.size() != n_params)
            {
                logError() << "line " << line_nb << " is misformatted";
                logError() << "  expecting" << n_params << "parameters (" << line_split.size() << "found)";
                logFatal(__FILE__, __func__, __LINE__);
            }
            else
            {
                for (size_t idx_param = 0; idx_param < n_params; idx_param++)
                {
           
                    bool conversion_success = true;
                    double val = swString::stringToDouble(line_split[idx_param], conversion_success);
                    if (conversion_success == false)
                    {
                        logError() << "line is misformatted: format conversion";
                        swIO::errorReadingLine(file_r, file_name, line_txt, line_nb, __LINE__, __func__, __FILE__);
                    }
                    param_array[idx_param].push_back(val);
                }
            }
            ++line_nb;
        }
        //====== Variable case

        // Warning if there are additional records in file
        if ((!file_r.eof()) && (swIO::readLine(file_r) != ""))
        {
            logWarn() << "Additional lines found in file " << file_name << " will be ignored.";
        }

        // Close file
        file_r.close();
    }
    logInfo() << "...done.";
    //==== read file
}


double swIO::timeFromFileName(std::string filename)
/**
 * @brief Compute time (in seconds) from the name of a file
 * @param[in] file_name : name of the file 
 * @return time in seconds (double)
 */
{
    std::string time_string;
    if (filename.substr(0, 12).compare("ResultHydro_") == 0)
    {
        time_string = filename.substr(filename.find("ResultHydro_") + 13, 17);
    }
    else if (filename.substr(0, 11).compare("ResultTime_") == 0)
    {
        time_string = filename.substr(filename.find("ResultTime_") + 12, 17);
    }
    else
    {
        logError() << "File format not recognized :" << filename;
        logError() << "Known formats are : ResultHydro_*, ResultTime_*.txt";
        logFatal(__FILE__, __func__, __LINE__);
    }

    bool conversion_success = true;
    size_t days = swString::stringToUInt(time_string.substr(0,4), conversion_success);
    size_t hours = swString::stringToUInt(time_string.substr(5,2), conversion_success);
    size_t minutes = swString::stringToUInt(time_string.substr(8,2), conversion_success);
    size_t seconds = swString::stringToUInt(time_string.substr(11,2), conversion_success);
    size_t milliseconds = swString::stringToUInt(time_string.substr(14,3), conversion_success);

    return (milliseconds / 1000. +
            seconds + 60 * (minutes + 60 * (hours + 24 * days)));
}


std::string swIO::readLine(std::ifstream &file)
/**
 * @brief read the next line of an open stream
 * @param[in] file
 * @return std::string content of the line
 */
{
    std::string line_txt;
    std::getline(file, line_txt);
    line_txt = std::regex_replace(line_txt, std::regex("(\\r|\\n)+"), "");
    return line_txt;
}


std::vector<std::string> swIO::stringSplit(std::string const& str) 
/**
 * @brief Splits a string into several words according to predefined separators
 * @param[in] str : the string to be split
 * @return a vector of strings 
*/
{
    std::vector<std::string> tokens = {};
    
    //-- Define characters that will be understood as being keywords"  
    std::regex rgx("[a-zA-Z0-9-_.+]+");
    
    auto begin = std::sregex_iterator(str.begin(), str.end(), rgx);
    auto end = std::sregex_iterator();
    
    for (auto it = begin; it != end; ++it) 
	{
        tokens.push_back(it->str());
    }
    
    return tokens;
}


bool swIO::writeTxtLine(std::ofstream& output_stream, const std::string file_line)
/**
 * @brief write line in an (open) streamfile
 * @param[out] output_stream: stream in which le line is written
 * @param[in] fileline: string to we written
 */
{
    output_stream << file_line;
    output_stream.flush();
    return EXIT_SUCCESS;
}
