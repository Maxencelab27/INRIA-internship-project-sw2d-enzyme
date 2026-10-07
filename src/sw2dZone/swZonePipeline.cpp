#include "swZonePipeline.h"

#include <filesystem>

#include "../common/geometry/swCell.h"
#include "../common/geometry/swEdgeBoundary.h"
#include "../common/geometry/swEdgeIntern.h"
#include "../common/geometry/swMesh.h"
#include "../common/swFile.h"
#include "../common/swIO.h"
#include "../common/swLog.h"
#include "swZone.h"

swZonePipeline::swZonePipeline(swMesh *mesh, swIO *io)
/**
 * @brief Construct a new swZonePipeline::swZonePipeline object
 * 
 * @param[in] mesh pointer to swMesh object
 * @param[in] io pointer to swIO object
 */
{
    m_mesh = mesh;
    m_io = io;
}

bool swZonePipeline::assignElements()
/**
 * @brief Assign cells, boundary edges and internal edges to their corresponding zone
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    for (swZone *zone : list_zone)
    {
        if (zone->assignElements() != EXIT_SUCCESS)
        {
            logError() << "Error in zone->assignElements() for zone" << zone->getName();
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

bool swZonePipeline::computeParam()
/**
 * @brief Check variables consistency and write the output files
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    //-- check consistency
    bool incoherence = false;
    for (size_t i = 0; i < var_per_file.size(); i++)
    {
        std::vector<std::string> variables = var_per_file[i];
        for (size_t j = 0; j < variables.size(); j++)
        {
            bool defined = false;
            for (std::string &var : default_variables_names)
            {
                if (var == variables[j])
                {
                    defined = true;
                    break;
                }
            }

            if (defined == false)
            {
                logError() << "Default value (Block2) have to be defined for variable" << variables[j];
                incoherence = true;
            }
        }
    }

    if (incoherence == true)
    {
        logError() << "parameters are incoherent";
        return EXIT_FAILURE;
    }
    //== check consistency

    //-- write files
    for (size_t k = 0; k < generated_files.size(); k++)
    {
        if (exportValue(k) != EXIT_SUCCESS)
        {
            logError() << "Error while writing" << generated_files[k];
            return EXIT_FAILURE;
        }
    }
    //== write files

    return EXIT_SUCCESS;
}

bool swZonePipeline::exportValue(size_t file_idx)
/**
 * @brief Create the values table for one file to write and call writeFile to write the file
 * 
 * @param[in] file_idx index of the file in generated_file attribute
 * @return bool EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    std::string file_name = swFile::input_folder + generated_files[file_idx];
    std::vector<std::string> file_vars = var_per_file[file_idx];

    size_t nb_elem = 0;

    if(file_elements.find(file_name) == file_elements.end())
    {
        logError() << "Unknown file :" << file_name;
        logFatal(__FILE__, __func__, __LINE__);
    }
    std::string elem_type = file_elements[file_name];

    if (elem_type == "cell")
    {
        nb_elem = m_mesh->getCellsList().size();
    }
    else if (elem_type == "edge_bnd")
    {
        nb_elem = m_mesh->getEdgeBoundaryList().size();
    }
    else if (elem_type == "edge_int")
    {
        nb_elem = m_mesh->getEdgeInternList().size();
    }

    size_t nb_col = file_vars.size();

    std::vector<std::vector<double>> values_table(nb_elem, std::vector<double>(nb_col));

    //-- set default values
    std::vector<double> file_default_values(nb_col);
    for (size_t col = 0; col < nb_col; ++col)
    {
        std::string var = file_vars[col];
        for (size_t idx_in_default = 0; idx_in_default < default_variables_names.size(); idx_in_default++)
        {
            if (default_variables_names[idx_in_default] == var)
            {
                file_default_values[col] = default_variables_values[idx_in_default];
            }
        }
        for (size_t line = 0; line < nb_elem; ++line)
        {
            for (size_t k = 0; k < default_variables_names.size(); k++)
            {
                if (default_variables_names[k] == var)
                {
                    values_table[line][col] = default_variables_values[k];
                    break;
                }
            }
        }
    }
    //== set default values

    for (swZone *zone : list_zone)
    {
        std::vector<std::string> zone_vars = zone->getVariables();
        for (size_t var_idx_in_zone = 0; var_idx_in_zone < zone_vars.size(); var_idx_in_zone++)
        {
            std::string zone_var = zone_vars[var_idx_in_zone];
            for (size_t var_idx_in_file = 0; var_idx_in_file < file_vars.size(); var_idx_in_file++)
            {
                if (file_vars[var_idx_in_file] == zone_var)
                {
                    if (elem_type == "cell")
                    {
                        std::vector<swCell *> zone_elems = zone->getCellList();
                        for (swCell *zone_elem : zone_elems)
                        {
                            values_table[zone_elem->getIdx()][var_idx_in_file] = zone->getValue(var_idx_in_zone);
                        }
                    }
                    else if (elem_type == "edge_bnd")
                    {
                        std::vector<swEdgeBoundary *> zone_elems = zone->getEdgeBoundaryList();
                        for (swEdgeBoundary *zone_elem : zone_elems)
                        {
                            values_table[zone_elem->getIdx()][var_idx_in_file] = zone->getValue(var_idx_in_zone);
                        }
                    }
                    else if (elem_type == "edge_int")
                    {
                        std::vector<swEdgeIntern *> zone_elems = zone->getEdgeInternList();
                        for (swEdgeIntern *zone_elem : zone_elems)
                        {
                            values_table[zone_elem->getIdx()][var_idx_in_file] = zone->getValue(var_idx_in_zone);
                        }
                    }
                }
            }
        }
    }

    if (writeFile(file_name, values_table, file_vars, file_default_values) != EXIT_SUCCESS)
    {
        logError() << "Error in writeFile";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

bool swZonePipeline::readInput()
/**
 * @brief Read the zone input file to set attributes values
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    std::string filename = swFile::input_zone_file;
    std::ifstream file(m_io->absolutePath(filename));
    if (file.is_open() == false)
    {
        logError() << "    Unable to open" << filename;
        logError() << __FILE__ << __func__ << __LINE__;
        return EXIT_FAILURE;
    }

    logInfo() << "Importing parameters for" << filename;

    std::size_t line_nb = 0;
    bool block1_defined = false;
    bool block2_defined = false;
    bool block3_defined = false;
    bool block4_defined = false;

    while (!file.eof())
    {
        std::string line_txt = swIO::readLine(file);
        ++line_nb;

        //-- Block1: zones list BEGIN (compulsory line)
        if (line_txt == "Block1: zones list BEGIN (compulsory line)")
        {
            line_txt = swIO::readLine(file);
            ++line_nb;
            while (!(line_txt == "Block1: zones list END (compulsory line)"))
            {
                if (!(line_txt.empty() || (line_txt[0] == '#') || (line_txt[0]== '=')))
                {
                    swZone *zone = swZone::create(m_mesh, line_txt);
                    if (zone == nullptr)
                    {
                        logError() << "Incoherent data to construct a zone (line" << line_nb <<")";
                        logError() << __FILE__ << __func__ << __LINE__;
                        return EXIT_FAILURE;
                    }
                    list_zone.push_back(zone);
                }

                line_txt = swIO::readLine(file);
                ++line_nb;
            }
            block1_defined = true;
        }
        //== Block1: zones list END (compulsory line)

        //-- Block2: default parameter value BEGIN (compulsory line)
        if (line_txt == "Block2: default parameter value BEGIN (compulsory line)")
        {
            line_txt = swIO::readLine(file);
            ++line_nb;
            while (!(line_txt == "Block2: default parameter value END (compulsory line)"))
            {
                if (!(line_txt.empty() || (line_txt[0] == '#') || (line_txt[0]== '=')))
                {
                    const std::vector<std::string>& line_split = swString::splitStr(line_txt, swString::regexWhiteSpace());
                    bool conversion_success = true;
                    double value = swString::stringToDouble(line_split[1], conversion_success);
                    if (conversion_success == false)
                    {
                        logError() << "line is misformatted (line" << line_nb <<")";
                        logError() << __FILE__ << __func__ << __LINE__;
                        return EXIT_FAILURE;
                    }

                    default_variables_names.push_back(line_split[0]);
                    default_variables_values.push_back(value);
                }

                line_txt = swIO::readLine(file);
                ++line_nb;
            }
            block2_defined = true;
        }
        //== Block2: default parameter value END (compulsory line)

        //-- Block3: linking variable-zone BEGIN (compulsory line)
        if (line_txt == "Block3: linking variable-zone BEGIN (compulsory line)")
        {
            line_txt = swIO::readLine(file);
            ++line_nb;
            while (!(line_txt == "Block3: linking variable-zone END (compulsory line)"))
            {
                if (!(line_txt.empty() || (line_txt[0] == '#') || (line_txt[0]== '=')))
                {
                    const std::vector<std::string>& line_split = swString::splitStr(line_txt, swString::regexWhiteSpace());
                    std::string var = line_split[0];
                    std::string zone_name = line_split[2];

                    bool conversion_success = true;
                    double value = swString::stringToDouble(line_split[1], conversion_success);
                    if (conversion_success == false)
                    {
                        logError() << "line is misformatted (line" << line_nb <<")";
                        logError() << __FILE__ << __func__ << __LINE__;
                        return EXIT_FAILURE;
                    }

                    int id = -1;
                    for (size_t i = 0; i < list_zone.size(); ++i)
                    {
                        if (zone_name == list_zone[i]->getName())
                        {
                            id = i;
                        }
                    }
                    if (id == -1)
                    {
                        logError() << "Line" << line_nb;
                        logError() << "Impossible to define parameter for zone" << zone_name;
                        logError() << "This zone has not been defined in Block1";
                        return EXIT_FAILURE;
                    }

                    if (list_zone[id]->addVariable(var, value, default_variables_names) != EXIT_SUCCESS)
                    {
                        logError() << "Error in zone->addVariable()";
                        logError() << var << value << default_variables_names;
                        logError() << __FILE__ << __func__ << __LINE__;
                        return EXIT_FAILURE;
                    }
                }

                line_txt = swIO::readLine(file);
                ++line_nb;
            }
            block3_defined = true;
        }
        //== Block3: linking variable-zone END (compulsory line)

        //-- Block4: output file name BEGIN (compulsory line)
        if (line_txt == "Block4: output file name BEGIN (compulsory line)")
        {
            line_txt = swIO::readLine(file);
            ++line_nb;
            while (!(line_txt == "Block4: output file name END (compulsory line)"))
            {
                if (!(line_txt.empty() || (line_txt[0] == '#') || (line_txt[0]== '=')))
                {
                    const std::vector<std::string>& line_split = swString::splitStr(line_txt, swString::regexWhiteSpace());

                    generated_files.push_back(line_split[0]);
                    std::vector<std::string> variables;
                    for (size_t k=1; k<line_split.size(); k++)
                    {
                        variables.push_back(line_split[k]);
                    }
                    var_per_file.push_back(variables);
                }

                line_txt = swIO::readLine(file);
                ++line_nb;
            }
            block4_defined = true;
        }
        //== Block4: output file name END (compulsory line)
    }
    file.close();

    ///-- Check blocks have all been read
    if (!(block1_defined && block2_defined && block4_defined))
    {
        logError() << "Missing zone blocks in" << swFile::input_zone_file;
        logError() << "Blocks read (1, 2, 3, 4): " << block1_defined << block2_defined << block3_defined << block4_defined;
        return EXIT_FAILURE; 
    }
    else if (!block3_defined)
    {
        logWarn() << "Warning, block 3 not defined (only uniform values considered)";
        return EXIT_SUCCESS;
    }
    else 
    {    
        logInfo() << "  Number of zones:" << list_zone.size();
        logInfo() << "...Done";
        return EXIT_SUCCESS;
    }
    //== Check blocks have all been read
}

bool swZonePipeline::writeFile(std::string filename, std::vector<std::vector<double>> &values, std::vector<std::string> &variable_names, std::vector<double> &default_values_list)
/**
 * @brief Write a zone output file
 * 
 * @param[in] filename file name
 * @param[in] values tab of values to write in the file
 * @param[in] variable_names list of the variable names
 * @param[in] default_values_list list of the variable default values
 * @return EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    logInfo() << "Writing" << filename;

    std::ofstream file;
    file.open(m_io->absolutePath(filename));
    if (file.is_open() == false)
    {
        logError() << "    Unable to open" << filename;
        logError() << __FILE__ << __func__ << __LINE__;
        return EXIT_FAILURE;
    }

    file << "Unif\t0" << std::endl;
    file << "==== Default Param" << std::endl;
    file << variable_names[0];
    for (size_t col = 1; col < variable_names.size(); ++col)
    {
        file << "\t" << variable_names[col];
    }
    file << std::endl;
    file << default_values_list[0];
    for (size_t col = 1; col < default_values_list.size(); ++col)
    {
        file << "\t" << default_values_list[col];
    }
    file << std::endl;
    file << "==== Distrib" << std::endl;

    for (size_t line = 0; line < values.size(); ++line)
    {
        file << values[line][0];
        for (size_t col = 1; col < values[line].size(); ++col)
        {
            file << "\t" << values[line][col];
        }
        file << std::endl;
    }

    file.flush();
    file.close();

    logInfo() << "... done";
    return EXIT_SUCCESS;
}
