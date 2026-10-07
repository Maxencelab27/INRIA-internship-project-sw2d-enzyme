#include "swPostProcessPipeline.h"

#include <cmath>
#include <filesystem>

#include "../common/swFile.h"
#include "../common/swIO.h"
#include "../common/swLog.h"

swPostProcessPipeline::swPostProcessPipeline(swIO *io, swPostProcessParam *param)
/**
 * @brief Construct a new swPostProcessPipeline::swPostProcessPipeline object
 * 
 * @param[in] io pointer to the swIO object
 * @param[in] param pointer to the swPostProcessParam object
 */
{
    m_io = io;
    m_param = param;
}

bool swPostProcessPipeline::computeMax()
/**
 * @brief Read simulation results and compute the maximum value and its date for each field and cell
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    //-- Import results
    x_cells.clear();
    y_cells.clear();
    fields.clear();
    times.clear();
    values.clear();

    importResults(m_param->getInputFormat(), m_param->getInputMax());
    //== Import results

    //-- Initialize max and dates of max tables
    size_t nb_cells = x_cells.size();
    size_t nb_fields = fields.size();
    std::vector<std::vector<double>> max_table = values[0];
    std::vector<std::vector<double>> dates_of_max_table(nb_fields, std::vector<double>(nb_cells, times[0]));
    //== Initialize max and dates of max tables

    //-- Compute max
    logInfo() << "Computing max and date of max";
    for (size_t c = 0; c < nb_cells; c++)
    {
        for (size_t f = 0; f < nb_fields; f++)
        {
            std::string field = fields[f];

            for (size_t t = 0; t < times.size(); t++)
            {
                if (field.find("_norm") != std::string::npos) // vector field
                {
                    if (values[t][f][c] > max_table[f][c])
                    {
                        max_table[f-2][c] = values[t][f-2][c]; // x component
                        max_table[f-1][c] = values[t][f-1][c]; // y component
                        max_table[f][c]   = values[t][f][c];   // norm

                        dates_of_max_table[f-2][c] = times[t]; // x component
                        dates_of_max_table[f-1][c] = times[t]; // y component
                        dates_of_max_table[f][c]   = times[t]; // norm
                    }
                }
                else if ((field.find("_x") == std::string::npos) &&
                            (field.find("_y") == std::string::npos)) // scalar field
                {
                    if (values[t][f][c] > max_table[f][c])
                    {
                        max_table[f][c] = values[t][f][c];
                        dates_of_max_table[f][c] = times[t];
                    }
                }
            }
        }
    }
    logInfo() << "...done";
    //== Compute max

    //-- Write output files
    // add x, y and zb fields
    fields.insert(fields.begin(), "bottom_elevation");
    fields.insert(fields.begin(), "y");
    fields.insert(fields.begin(), "x");
    max_table.insert(max_table.begin(), zb_cells);
    max_table.insert(max_table.begin(), y_cells);
    max_table.insert(max_table.begin(), x_cells);
    dates_of_max_table.insert(dates_of_max_table.begin(), zb_cells);
    dates_of_max_table.insert(dates_of_max_table.begin(), y_cells);
    dates_of_max_table.insert(dates_of_max_table.begin(), x_cells);

    // write
    if (swIO::createTxtFileMapDecimal(m_io->absolutePath(swFile::output_map_max_file), fields, max_table, 10) != EXIT_SUCCESS)
    {
        logError() << "Failure during writting max file" << swFile::output_map_max_file;
        return EXIT_FAILURE;
    }
    if (swIO::createTxtFileMapDecimal(m_io->absolutePath(swFile::output_map_max_date_file), fields, dates_of_max_table, 10) != EXIT_SUCCESS)
    {
        logError() << "Failure during writting dates of max file" << swFile::output_map_max_file;
        return EXIT_FAILURE;
    }
    //== Write output files

    return EXIT_SUCCESS;
}

bool swPostProcessPipeline::computeSubmersions()
/**
 * @brief Read simulation results and compute the date of first submersion and submersion durations for each cell
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    //-- Import results
    x_cells.clear();
    y_cells.clear();
    fields.clear();
    times.clear();
    values.clear();

    importResults(m_param->getInputFormat(), m_param->getInputSubmersion());
    //== Import results

    //-- Get water depth values
    size_t i_h = 0;
    while (fields[i_h] != "water_depth")
    {
        i_h++;
    }

    std::vector<std::vector<double>> water_depths;
    for (size_t k = 0; k < values.size(); k++)
    {
        water_depths.push_back(values[k][i_h]);
    }
    //== Get water depth values

    //-- Sort values
    size_t i_min;
    for (size_t i = 0; i < times.size(); i++)
    {
        i_min = i;
        for (size_t j = i; j < times.size(); j++)
        {
            if (times[j] < times[i_min])
            {
                i_min = j;
            }
        }
        std::swap(times[i], times[i_min]);
        std::swap(water_depths[i], water_depths[i_min]);
    }
    //== Sort values

    //-- Compute submersion durations
    logInfo() << "Computing submersion durations";
    double h_min = m_param->getHmin();
    size_t nb_cells = x_cells.size();

    std::vector<double> first_sub(nb_cells, -999.);
    std::vector<double> last_sub(nb_cells, -999.);
    std::vector<double> last_drying(nb_cells, -999.);
    std::vector<double> max_sub(nb_cells, 0.);
    std::vector<double> total_sub(nb_cells, 0.);


    for (size_t n = 0; n < nb_cells; n++)
    {
        if (water_depths[0][n] > h_min) // if cell initially submerged
        {
            first_sub[n] = times[0];
            last_sub[n] = times[0];
        }
    }

    //compute max submersion and first flooding
    for (size_t t = 1; t < times.size(); t++)
    {
        for (size_t n = 0; n < nb_cells; n++)
        {
            if ((water_depths[t-1][n] <= h_min) && (water_depths[t][n] > h_min)) // if cell get submerged
            {
                last_sub[n] = times[t];
                if (first_sub[n] == -999.) // if it's the first time
                {
                    first_sub[n] = times[t];
                }
            }
            if ((water_depths[t-1][n] > h_min) && (water_depths[t][n] <= h_min)) // if cell get dry
            {
                last_drying[n] = times[t];
                double sub_duration = last_drying[n] - last_sub[n];
                if (sub_duration > max_sub[n]) // check if it's the max submersion duration
                {
                    max_sub[n] = sub_duration;
                }
            }
        }
    }
    for (size_t n = 0; n < nb_cells; n++)
    {
        if (water_depths[times.size()-1][n] > h_min) //if cell is submerged at the end
        {
            last_drying[n] = times[times.size()-1];
            double sub_duration = last_drying[n] - last_sub[n];
            if (sub_duration > max_sub[n]) // check if it's the max submersion duration
            {
                max_sub[n] = sub_duration;
            }
        }
    }

    // compute total submersion
    for (size_t n = 0; n < nb_cells; n++)
    {
        if (first_sub[n] != 999.)
        {
            total_sub[n] = last_drying[n] - first_sub[n];
        }
    }
    logInfo() << "...done";
    //== Compute submersion durations

    //-- Write output files
    std::vector<std::string> var_names = {"x", "y", "zb", "first_submersion_date", "max_submersion_duration", "total_submersion_duration"};
    std::vector<std::vector<double>> var_values = {x_cells, y_cells, zb_cells, first_sub, max_sub, total_sub};

    if (swIO::createTxtFileMapDecimal(m_io->absolutePath(swFile::output_map_submersion_durations), var_names, var_values, 10) != EXIT_SUCCESS)
    {
        logError() << "Failure during writting submersion durations file" << swFile::output_map_submersion_durations;
        return EXIT_FAILURE;
    }
    //== Write output files

    return EXIT_SUCCESS;
}

bool swPostProcessPipeline::importBinResults()
/**
 * @brief read the simulation results from a binary file and set the attributes values
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    // Open file
    std::string filename = swFile::output_map_binary_file;
    logInfo() << "Reading file:" << filename;

    std::ifstream file(m_io->absolutePath(filename), std::ios::binary);
    if (!file.is_open())
    {
        logError() << "  Unable to open file" << filename;
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    //== Open file

    //-- Read geometry
    // Read the first line (version number)
    size_t version_lenght;
    file.read(reinterpret_cast<char*>(&version_lenght), sizeof(version_lenght));
    std::string read_version(version_lenght, '\0');
    file.read(&read_version[0], version_lenght);

    size_t nb_cell;
    file.read(reinterpret_cast<char*>(&nb_cell), sizeof(nb_cell));
    size_t nb_edges_boundary;
    file.read(reinterpret_cast<char*>(&nb_edges_boundary), sizeof(nb_edges_boundary));
    size_t nb_edges_intern;
    file.read(reinterpret_cast<char*>(&nb_edges_intern), sizeof(nb_edges_intern));
    size_t nb_node;
    file.read(reinterpret_cast<char*>(&nb_node), sizeof(nb_node));
    size_t nb_boundaries;
    file.read(reinterpret_cast<char*>(&nb_boundaries), sizeof(nb_boundaries));

    // Read cells
    std::vector<size_t> cell_nb_node;
    std::vector<size_t> cell_nb_edge_bnd;
    std::vector<size_t> cell_nb_edge_int;
    std::vector<size_t> cell_nb_neigh;

    for (size_t n = 0; n < nb_cell; ++n)
    {
        size_t nb_node;
        size_t nb_edge_bnd;
        size_t nb_edge_int;
        size_t nb_neigh;
        double xg;
        double yg;
        double zb;
        file.seekg(sizeof(size_t) * 2, std::ios_base::cur); //skip unused data
        file.read(reinterpret_cast<char*>(&nb_node), sizeof(nb_node));
        file.read(reinterpret_cast<char*>(&nb_edge_bnd), sizeof(nb_edge_bnd));
        file.read(reinterpret_cast<char*>(&nb_edge_int), sizeof(nb_edge_int));
        file.read(reinterpret_cast<char*>(&nb_neigh), sizeof(nb_neigh));
        file.seekg(sizeof(double) * 2, std::ios_base::cur); //skip unused data
        file.read(reinterpret_cast<char*>(&xg), sizeof(xg));
        file.read(reinterpret_cast<char*>(&yg), sizeof(yg));
        file.read(reinterpret_cast<char*>(&zb), sizeof(zb));

        x_cells.push_back(xg);
        y_cells.push_back(yg);
        zb_cells.push_back(zb);
        cell_nb_node.push_back(nb_node);
        cell_nb_edge_bnd.push_back(nb_edge_bnd);
        cell_nb_edge_int.push_back(nb_edge_int);
        cell_nb_neigh.push_back(nb_neigh);
    }

    // Read nodes
    std::vector<size_t> node_nb_cell;
    for (size_t p = 0; p < nb_node; ++p)
    {
        size_t ncelp;
        file.seekg(sizeof(size_t) * 2 + sizeof(double) * 3, std::ios_base::cur); //skip unused data
        file.read(reinterpret_cast<char*>(&ncelp), sizeof(ncelp));
        file.seekg(sizeof(size_t), std::ios_base::cur); //skip unused data
        node_nb_cell.push_back(ncelp);
    }

    // Skip boundary edges
    file.seekg((sizeof(size_t)*4 + sizeof(int) + sizeof(double)*7) * nb_edges_boundary, std::ios_base::cur);

    // Skip internal edges
    file.seekg((sizeof(size_t)*5 + sizeof(double)*7) * nb_edges_intern, std::ios_base::cur);

    // Skip boundaries
    for (size_t i_bnd = 0; i_bnd < nb_boundaries; i_bnd++)
    {
        size_t bnd_size;
        file.read(reinterpret_cast<char*>(&bnd_size), sizeof(bnd_size));
        file.seekg(sizeof(size_t) * bnd_size, std::ios_base::cur);
    }
    file.seekg(sizeof(int) * nb_boundaries, std::ios_base::cur);

    // Skip cell connectivity
    for (size_t k = 0; k < nb_cell; k++)
    {
        file.seekg(sizeof(size_t) * (cell_nb_node[k] + cell_nb_edge_bnd[k] + cell_nb_edge_int[k] + cell_nb_neigh[k]), std::ios_base::cur);
    }

    // Skip node connectivity
    for (size_t k = 0; k < nb_node; k++)
    {
        file.seekg(sizeof(size_t) * node_nb_cell[k], std::ios_base::cur);
    }
    //== Read geometry

    //-- Read variable names
    size_t nb_var;
    file.read(reinterpret_cast<char*>(&nb_var), sizeof(nb_var));
    
    for (size_t k = 0; k < nb_var; k++)
    {
        size_t name_length;
        file.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
        std::string name(name_length, '\0');
        file.read(&name[0], name_length);
        if (k != 0) // skip zb
        {
            fields.push_back(name);
        }
    }
    //== Read variable names

    //-- Read time information
    double t_0;
    double t_max;
    double dt_map;
    file.read(reinterpret_cast<char*>(&t_0), sizeof(t_0));
    file.read(reinterpret_cast<char*>(&t_max), sizeof(t_max));
    file.read(reinterpret_cast<char*>(&dt_map), sizeof(dt_map));
    //== Read time information

    //-- Read results
    size_t nb_fields = fields.size();
    double t = t_0;
    
    while (t <= t_max)
    {
        file.read(reinterpret_cast<char*>(&t), sizeof(t));
        times.push_back(t);

        std::vector<std::vector<double>> time_values(nb_fields, std::vector<double>(nb_cell, -999.));
        for (size_t c = 0; c < nb_cell; c++)
        {
            for (size_t f = 0; f < nb_fields + 1; f++)
            {
                double value;
                file.read(reinterpret_cast<char*>(&value), sizeof(value));
                if (f != 0) // skip zb
                {
                    time_values[f-1][c] = value;
                }
            }
        }
        values.push_back(time_values);

        if (t == t_max) {break;}
    }
    //== Read results

    return EXIT_SUCCESS;
}

bool swPostProcessPipeline::importResults(std::string format_option, std::string files_option)
/**
 * @brief read the simulation results set the attributes values
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    if (format_option == "txt")
    {
        if (((files_option == "all") && !result_hydro_imported && result_time_imported) ||
            ((files_option == "resulthydro") && !result_hydro_imported))
        {
            if (importTxtResults("resulthydro") != EXIT_SUCCESS)
            {
                logError() << "Failure while import text result files";
                return EXIT_FAILURE;
            }
            result_hydro_imported = true;
        }
        if (((files_option == "all") && result_hydro_imported && !result_time_imported) ||
            ((files_option == "resulttime") && !result_time_imported))
        {
            if (importTxtResults("resulttime") != EXIT_SUCCESS)
            {
                logError() << "Failure while import text result files";
                return EXIT_FAILURE;
            }
            result_time_imported = true;
        }
        if ((files_option == "all") && !result_hydro_imported && !result_time_imported)
        {
        {
            if (importTxtResults("all") != EXIT_SUCCESS)
            {
                logError() << "Failure while import text result files";
                return EXIT_FAILURE;
            }
            result_hydro_imported = true;
            result_time_imported = true;
        }
        }
    }
    if ((format_option == "bin") && !result_bin_imported)
    {
        if (importBinResults() != EXIT_SUCCESS)
        {
            logError() << "Failure while import binary result files";
            return EXIT_FAILURE;
        }
        result_bin_imported = true;
    }

    return EXIT_SUCCESS;
}

bool swPostProcessPipeline::importTxtResults(std::string files_option)
/**
 * @brief read the simulation results from a text files and set the attributes values
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    //-- Get the list of file names
    logInfo() << "Scan output folder:" << swFile::output_folder;

    std::vector<std::string> files_list;
    if (std::filesystem::exists(m_io->absolutePath(swFile::output_folder)) == false)
    {
        logError() << "Missing folder:" << swFile::output_folder;
        logFatal(__FILE__, __func__, __LINE__);
    }
    for (const auto& entry : std::filesystem::directory_iterator(m_io->absolutePath(swFile::output_folder)))
    {
        std::string filepath = entry.path().filename().u8string();
        if ((filepath.substr(0, 13).compare("ResultHydro_p") == 0) &&
            ((files_option == "resulthydro") || (files_option == "all")))
        {
            files_list.push_back(filepath);
        }
        if ((filepath.substr(0, 12).compare("ResultTime_p") == 0) &&
            ((files_option == "resulttime") || (files_option == "all")))
        {
            files_list.push_back(filepath);
        }
    }
    if (files_list.size() == 0)
    {
        logError() << "No text result file found in output folder :" << swFile::output_folder;
        logFatal(__FILE__, __func__, __LINE__);
    }

    logInfo() << "...done";
    //== Get the list of file names

    logInfo() << "Reading files from:" << swFile::output_folder;

    //-- Read first file (for x, y, field names and nb_cells)
    std::string filename = files_list[0];

    // open file
    std::string filepath =  m_io->absolutePath(swFile::output_folder + filename);
    std::ifstream file(filepath);

    // import fields
    std::string line_txt = swIO::readLine(file);
    std::vector<std::string> line_split = swString::splitStr(line_txt, swString::regexWhiteSpace());
    for (size_t f = 3; f < line_split.size(); f++)
    {
        fields.push_back(line_split[f]);
    }

    line_txt = swIO::readLine(file);
    while (!file.eof())
    {
        std::vector<std::string> cell_values_txt = swString::splitStr(line_txt, swString::regexWhiteSpace());
        bool conversion_success = true;

        //import x_cells, y_cells and zb_cells
        double x = swString::stringToDouble(cell_values_txt[0], conversion_success);
        double y = swString::stringToDouble(cell_values_txt[1], conversion_success);
        double zb = swString::stringToDouble(cell_values_txt[2], conversion_success);
        x_cells.push_back(x);
        y_cells.push_back(y);
        zb_cells.push_back(zb);

        line_txt = swIO::readLine(file);
    }

    size_t nb_fields = fields.size();
    size_t nb_cells = x_cells.size();
    //== Read first file (for x, y, field names and nb_cells)

    //-- Read files
    for (std::string filename : files_list)
    {
        // import time
        times.push_back(swIO::timeFromFileName(filename));
        
        // open file
        std::string filepath =  m_io->absolutePath(swFile::output_folder + filename);
        std::ifstream file(filepath);

        // read header line
        line_txt = swIO::readLine(file);

        line_txt = swIO::readLine(file);
        std::vector<std::vector<double>> file_values(nb_fields, std::vector<double>(nb_cells, -999.));
        for (size_t c = 0; c < nb_cells; c++)
        {
            std::vector<std::string> cell_values_txt = swString::splitStr(line_txt, swString::regexWhiteSpace());
            bool conversion_success = true;

            // import values
            for (size_t f = 0; f < nb_fields; f++)
            {
                double val = swString::stringToDouble(cell_values_txt[f+3], conversion_success);
                file_values[f][c] = val;
            }

            line_txt = swIO::readLine(file);
        }
        values.push_back(file_values);

        logInfo() << "  file" << filename << "done";
    }
    logInfo() << "...done";
    //== Read files

    return EXIT_SUCCESS;
}