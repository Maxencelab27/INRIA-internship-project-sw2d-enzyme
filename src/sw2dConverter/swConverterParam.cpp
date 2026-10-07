#include "swConverterParam.h"

#include <filesystem>

#include "../common/swFile.h"
#include "../common/param.h"
#include "../common/swIO.h"
#include "../common/swLog.h"
#include "../common/swString.h"

swConverterParam::swConverterParam(swIO *io)
/**
 * @brief construct the converterParam class
 * 
 */
{
    m_io = io;
    logInfo() << "Simulation folder: " << m_io->project_folder << std::endl;
}


bool swConverterParam::checkConsistency()
/**
 * @brief check the consistency of the input parameters
 * @return 'true' if the input parameters are consistent
 */
{
    bool consistent = true;

    // consistensy for converter output file
    if ((std::filesystem::exists(m_io->absolutePath(swFile::input_geo_file)) == true) && (overwrite == false))
    {
        logError() << "The destination file exists";
        logError() << "  Specify keyword 'overwrite' to 1 (in the primary_files/input.conv file) if you really want to erase this file";
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    //

    // consistensy for default bc
    if (default_bc_defined == false)
    {
        logError() << "Incorrect inputs sw2dConverter";
        logError() << "  Option 'default_bc' has to be user defined";
        consistent = false;
    }

    // consistency for topography correction
    if (topo_correction == true)
    {
        if (topo_option == "")
        {
            logError() << "Topography correction is activated (user has either defined parameter 'topo_file' or 'topo_option'";
            logError() << "Topography file has been given but no topo_option parameter ";
            logError() << "  (mandatory, should be min, max, or mean)";
            consistent = false;
        }
    }
    //== check parameters consistency
    return consistent;
}


bool swConverterParam::correctTopography()
/**
 * @brief getter to the parameter 'topo_correction'
 * @return bool true is there is correction to be done, false if not
 */
{
    return topo_correction;
}


double swConverterParam::getAMin()
/**
 * @brief getter to the parameter 'a_min'
 */
{
    return a_min;
}


std::string swConverterParam::getBcFileName()
/**
 * @brief getter to the parameter 'bc_file'
 */
{
    return bc_file;
}


int swConverterParam::getDefaultBc()
/**
 * @brief getter to the parameter 'default_bc'
 */
{
    return default_bc;
}


std::string swConverterParam::getMeshFileName()
/**
 * @brief getter to the parameter 'mesh_file'
 */
{
    return mesh_file;
}


bool swConverterParam::getOverWrite()
/**
 * @brief getter to the parameter 'overwrite'
 */
{
    return overwrite;
}


std::string swConverterParam::getTopoFileName()
/**
 * @brief getter to the parameter 'topo_file'
 */
{
    return topo_file;
}


std::string swConverterParam::getTopoOption()
/**
 * @brief getter to the parameter 'topo_option'
 */
{
    return topo_option;
}


bool swConverterParam::getWriteLogs()
/**
 * @brief getter to the parameter 'write_logs'
 */
{
    return write_logs;
}


bool swConverterParam::readInputFile()
/**
 * @brief Read a converter input file and sets the corresponding converter parameters
 * 
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    // Open file
    std::ifstream file(m_io->absolutePath(swFile::input_converter_file));
    if (file.is_open() == false)
    {
        logError() << "  Unable to open:" << swFile::input_converter_file;
        logFatal(__FILE__, __func__, __LINE__);
    }

    size_t line_nb = 1;
    size_t exit_code = 0;
    while (!file.eof())
    {
        std::string line_txt = swIO::readLine(file);
        if (!(line_txt.empty() || (line_txt[0] == '#') || (line_txt[0]== '=')))
        {
            //-- Read line
            std::vector<std::string> line_split = swString::splitStr(line_txt, swString::regexWhiteSpace());
            if (line_split.size() != 2)
            {
                logError() << "  line" << line_nb << "is misformatted:" << line_txt;
                logError() << "    expected 2 elements : 1 parameter name and 1 parameter value";
                exit_code++;
            }
            //== Read
            else
            {
                //-- testing the allowed parameters
                std::string parameter_name = swString::toLower(line_split[0]);
                std::string parameter_value = line_split[1];
                if (parameter_name == "a_min")
                {
                    bool conversion_success = true;
                    a_min = swString::stringToDouble(parameter_value, conversion_success);
                    if (!((conversion_success == true) || (a_min > 0)))
                    {
                        logError() << "  unallowed value for parameter 'a_min'" << parameter_value;
                        logError() << "    'amin' should be a positive number";
                        exit_code++;
                    }
                }
                else if (parameter_name == "bc_file")
                {
                    bc_file = swFile::primary_folder + parameter_value;
                    if (std::filesystem::exists(m_io->absolutePath(bc_file)) == false)
                    {
                        logError() << "  unallowed value for parameter 'bc_file'" << parameter_value;
                        logError() << "    could not find file named" << bc_file;
                        exit_code++;
                    }
                }
                else if (parameter_name == "default_bc")
                {
                        bool conversion_success = true;
                        default_bc = swString::stringToInt(parameter_value, conversion_success);
                        default_bc_defined = true;
                        if ((conversion_success == false) || (default_bc >= 0))
                        {
                            logError() << "  unallowed value for parameter 'default_bc'" << parameter_value;
                            logError() << "    'default_bc' should be negative";
                            exit_code++;
                        }
                }
                else if (parameter_name == "mesh_file")
                {
                    mesh_file = swFile::primary_folder + parameter_value;
                    if (std::filesystem::exists(m_io->absolutePath(mesh_file)) == false)
                    {
                        logError() << "  unallowed value for parameter 'mesh_file'" << parameter_value;
                        logError() << "    could not find file named" << mesh_file;
                        exit_code++;
                    }
                }
                else if (parameter_name == "overwrite")
                {
                    bool conversion_success = true;
                    int val = swString::stringToInt(parameter_value, conversion_success);
                    if (!conversion_success || !((val == 0) || (val == 1)))
                    {
                        logError() << "  unallowed value for parameter 'overwrite'" << parameter_value;
                        logError() << "    allowed values for parameter 'overwrite'";
                        logError() << "    - 0\t do not overwrite existing files";
                        logError() << "    - 1\t please overwrite existing files";
                        exit_code++;
                    }
                    if (val == 0)
                    {
                        overwrite = false;
                    }
                    if (val == 1)
                    {
                        overwrite = true;
                    }
                }
                else if (parameter_name == "topo_file")
                {
                    std::string path = swFile::primary_folder + parameter_value;
                    if (std::filesystem::exists(m_io->absolutePath(path)) == false)
                    {
                        logError() << "  unallowed value for parameter 'topo_file'" << parameter_value;
                        logError() << "    could not find file named" << path;
                        exit_code++;
                    }
                    topo_file = path;
                    topo_correction = true;
                }
                else if (parameter_name == "topo_option")
                {
                    topo_option = swString::toLower(parameter_value);
                    if (!((topo_option == "min") || (topo_option == "max") || (topo_option == "mean")))
                    {
                        logError() << "  unallowed value for parameter 'topo_option'" << parameter_value;
                        logError() << "    allowed values for parameter 'topo_option'";
                        logError() << "    'min'\tno building exchange";
                        logError() << "    'max'\tcell_average building exchange";
                        logError() << "    'mean'\tcell_average building exchange";
                        exit_code++;
                    }
                    topo_correction = true;
                }
                else if (parameter_name == "write_logs")
                {
                    bool conversion_success = true;
                    int val = swString::stringToInt(parameter_value, conversion_success);
                    if (!conversion_success || !((val == 0) || (val == 1)))
                    {
                        logError() << "  unallowed value for parameter 'write_logs'" << parameter_value;
                        logError() << "    allowed values for parameter 'write_logs'";
                        logError() << "    - 0\t do not write log files";
                        logError() << "    - 1\t write log files";
                        exit_code++;
                    }
                    if (val == 0)
                    {
                        write_logs = false;
                    }
                    if (val == 1)
                    {
                        write_logs = true;
                    }
                }
                //== testing the allowed parameters
                else
                {
                    logError() << " Unknown parameter at line" << line_nb << ":" << line_split[0];
                    exit_code++;
                }
            }
        }

        ++line_nb;
    }

    if (exit_code > 0)
    {
        logError() << "  Unallowed value(s) in parameter file";
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


bool swConverterParam::writeLog()
/**
 * @brief Write all converter parameters into the log folder ()
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    // Open file
    std::ofstream file_parameters(m_io->absolutePath(swFile::log_input_conv_file));
    if (file_parameters.is_open() == false)
    {
        logError() << "    Can't create file " << swFile::log_input_conv_file;
        logFatal(__FILE__, __func__, __LINE__);
    }

    // Write parameters
    file_parameters << "a_min"         << "\t" << swString::numberToStr(a_min, 'e', 15) << std::endl;
    file_parameters << "bc_file"       << "\t" << bc_file << std::endl;
    file_parameters << "default_bc"    << "\t" << std::to_string(default_bc) << std::endl;
    file_parameters << "mesh_file"     << "\t" << mesh_file << std::endl;
    file_parameters << "overwrite"     << "\t" << std::to_string(overwrite) << std::endl;
    if (topo_correction == true)
    {
        file_parameters << "topo_file"     << "\t" << topo_file << std::endl;
        file_parameters << "topo_option"   << "\t" << topo_option << std::endl;
    }
    else
    {
        file_parameters << "topo_file"     << "\t" << "topography correction not activated = no file defined" << std::endl;
        file_parameters << "topo_option"   << "\t" << "topography correction not activated = no option defined" << std::endl;
    }

    file_parameters.close();

    return EXIT_SUCCESS;
}
