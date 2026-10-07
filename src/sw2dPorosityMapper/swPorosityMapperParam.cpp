#include "swPorosityMapperParam.h"

#include <filesystem>

#include "../common/swFile.h"
#include "../common/param.h"
#include "../common/swIO.h"
#include "../common/swLog.h"
#include "../common/swString.h"

swPorosityMapperParam::swPorosityMapperParam(swIO *io)
/**
 * @brief construct a swPorosityMapperParam object
 * 
 */
{
    m_io = io;
    logInfo() << "Simulation folder: " << m_io->project_folder << std::endl;
}


bool swPorosityMapperParam::checkConsistency()
/**
 * @brief check the consistency of the input parameters
 * @return 'true' if the input parameters are consistent
 */
{
    bool consistent = true;

    return consistent;
}

std::string swPorosityMapperParam::getInput()
/**
 * @brief return string of the input format for buildings
 */
{
    return input;
}

std::string swPorosityMapperParam::getModel()
/**
 * @brief return string of the porosity model
 */
{
    return model;
}

bool swPorosityMapperParam::readInputFile(std::string filepath)
/**
 * @brief Read a porosity mapper input file and sets the corresponding parameters
 * @param[in] filepath path to the porosity mapper input file
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    // Open file
    std::ifstream file(m_io->absolutePath(filepath));
    if (file.is_open() == false)
    {
        logError() << "  Unable to open:" << filepath;
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

                if (parameter_name == "model")
                {
                    model = swString::toLower(parameter_value);
                    if ((model != "ddp") &&
                        (model != "dip") &&
                        (model != "sp"))
                    {
                        logError() << "  Invalid value for parameter 'model_hyp'" << parameter_value;
                        logError() << "    allowed values for the 'model_hyp' parameter are";
                        logError() << "    - dip\tDepth Dependant Porosity model";
                        logError() << "    - dip\tDual Integral Porosity model";
                        logError() << "    - sp\tSingle Porosity model";
                        return EXIT_FAILURE;
                    }
                }
                else if (parameter_name == "input")
                {
                    input = swString::toLower(parameter_value);
                    if ((input != "footprint") &&
                        (input != "raster"))
                    {
                        logError() << "  Invalid value for parameter 'input'" << parameter_value;
                        logError() << "    allowed values for the 'input' parameter are";
                        logError() << "    - footprint\tInput file for buildings contains the building footprints";
                        logError() << "    - raster\tInput file for building is a raster";
                        return EXIT_FAILURE;
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


bool swPorosityMapperParam::writeLog(std::string filepath)
/**
 * @brief Write all porosity mapper parameters into the log folder
 * @param[in] filepath path to the log file for porosity mapper parameter
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    // Open file
    std::ofstream file_parameters(m_io->absolutePath(filepath));
    if (file_parameters.is_open() == false)
    {
        logError() << "    Can't create file " << filepath;
        logFatal(__FILE__, __func__, __LINE__);
    }

    // Write parameters
    file_parameters << "model"         << "\t" << model;
    file_parameters << "input"         << "\t" << input;

    file_parameters.close();

    return EXIT_SUCCESS;
}
