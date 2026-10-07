#include "swPostProcessParam.h"

#include <filesystem>

#include "../common/swFile.h"
#include "../common/param.h"
#include "../common/swIO.h"
#include "../common/swLog.h"
#include "../common/swString.h"

swPostProcessParam::swPostProcessParam(swIO *io)
/**
 * @brief construct a new postProcessParam::postProcessParam object
 */
{
    m_io = io;
    logInfo() << "Simulation folder: " << m_io->project_folder << std::endl;
}


bool swPostProcessParam::checkConsistency()
/**
 * @brief check the consistency of the input parameters
 * @return 'true' if the input parameters are consistent
 */
{
    bool consistent = true;

    return consistent;
}

bool swPostProcessParam::getComputeMax(void)
/**
 * @brief return the value of attribute 'compute_max'
 */
{
    return compute_max;
}

bool swPostProcessParam::getFirstSubmersion(void)
/**
 * @brief return the value of attribute 'first_submersion'
 */
{
    return first_submersion;
}

double swPostProcessParam::getHmin(void)
/**
 * @brief return the value of attribute 'h_min'
 */
{
    return h_min;
}

std::string swPostProcessParam::getInputFormat(void)
/**
 * @brief return the value of attribute 'input_format'
 */
{
    return input_format;
}

std::string swPostProcessParam::getInputMax(void)
/**
 * @brief return the value of attribute 'input_max'
 */
{
    return input_max;
}

std::string swPostProcessParam::getInputSubmersion(void)
/**
 * @brief return the value of attribute 'input_submersion'
 */
{
    return input_submersion;
}

std::string swPostProcessParam::getSubmersionDurations(void)
/**
 * @brief return the value of attribute 'submersion_durations'
 */
{
    return submersion_durations;
}

bool swPostProcessParam::readInputFile()
/**
 * @brief Read a postprocess input file and sets the corresponding postprocess parameters
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    // Open file
    std::ifstream file(m_io->absolutePath(swFile::input_postprocess_file));
    if (file.is_open() == false)
    {
        logError() << "  Unable to open:" << swFile::input_postprocess_file;
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
                if (parameter_name == "compute_max")
                {
                    exit_code += setcomputeMax(parameter_value);
                }
                else if (parameter_name == "first_submersion")
                {
                    exit_code += setFirstSubmersion(parameter_value);
                }
                else if (parameter_name == "h_min")
                {
                    exit_code += setHmin(parameter_value);
                }
                else if (parameter_name == "input_format")
                {
                    exit_code += setInputFormat(parameter_value);
                }
                else if (parameter_name == "input_max")
                {
                    exit_code += setInputMax(parameter_value);
                }
                else if (parameter_name == "input_submersion")
                {
                    exit_code += setInputSubmersion(parameter_value);
                }
                else if (parameter_name == "submersion_durations")
                {
                    exit_code += setSubmersionDurations(parameter_value);
                }
                //== testing the allowed parameters
                else
                {
                    logError() << "  - Line" << line_nb
                                << "Unknown parameter:" << line_split[0];
                }
            }
        }

        ++line_nb;
    }

    if (exit_code > 0)
    {
        logError() << "  Unallowed value(s) in parameter file";
        logFatal(__FILE__, __func__, __LINE__);
    }

    return EXIT_SUCCESS;
}

bool swPostProcessParam::setcomputeMax(const std::string& text)
/**
 * @brief set the value of attribute 'compute_max'
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    bool conversion_success = true; 
    int val = swString::stringToInt(text, conversion_success); 
    if (!conversion_success ||
        !((val == 0) || (val == 1)))
    {
        logError() << "  unallowed value for parameter 'compute_max'" << text;
        logError() << "    allowed values for parameter 'compute_max'";
        logError() << "    - 0\tno computation maximums";
        logError() << "    - 1\tcompute the maximums";
        return EXIT_FAILURE;
    }
    if (val == 0)
    {
        compute_max = false;
    }
    if (val == 1)
    {
        compute_max = true;
    }
    return EXIT_SUCCESS;
}

bool swPostProcessParam::setFirstSubmersion(const std::string& text)
/**
 * @brief set the value of attribute 'first_submersion'
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    bool conversion_success = true; 
    int val = swString::stringToInt(text, conversion_success); 
    if (!conversion_success ||
        !((val == 0) || (val == 1)))
    {
        logError() << "  unallowed value for parameter 'first_submersion'" << text;
        logError() << "    allowed values for parameter 'first_submersion'";
        logError() << "    - 0\tno computation of first submersion date";
        logError() << "    - 1\tcompute the date of first submersion";
        return EXIT_FAILURE;
    }
    if (val == 0)
    {
        first_submersion = false;
    }
    if (val == 1)
    {
        first_submersion = true;
    }
    return EXIT_SUCCESS;
}

bool swPostProcessParam::setHmin(const std::string& text)
/**
 * @brief set the value of attribute 'h_min'
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    bool conversion_success = true; 
    h_min = swString::stringToDouble(text, conversion_success);
    if (!conversion_success ||
        (h_min < 0))
    {
        logError() << "  unallowed value for parameter 'h_min'" << text;
        logError() << "    'h_min' should be positive";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

bool swPostProcessParam::setInputFormat(const std::string& text)
/**
 * @brief set the value of attribute 'input_format'
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    input_format = swString::toLower(text);
    if (!((input_format == "bin") ||
          (input_format == "txt")))
    {
        logError() << "  invalid value for parameter 'input_format'" << text;
        logError() << "  Allowed values are:";
        logError() << "    'bin'\t use binary output to compute max";
        logError() << "    'txt'\t use txt output to compute max";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

bool swPostProcessParam::setInputMax(const std::string& text)
/**
 * @brief set the value of attribute 'input_max'
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    input_max = swString::toLower(text);
    if (!((input_max == "resulthydro") ||
          (input_max == "resulttime") ||
          (input_max == "all")))
    {
        logError() << "  invalid value for parameter 'input_max'" << text;
        logError() << "  Allowed values are:";
        logError() << "    'all'\t Compute the maximum values for both ResultHydro**.txt and ResultTime**.txt files";
        logError() << "    'resulttime'\t  Compute the maximum values maps for ResultTime**.txt files only";
        logError() << "    'resulthydro'\t Compute the maximum values maps for ResultHydro**.txt files only";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

bool swPostProcessParam::setInputSubmersion(const std::string& text)
/**
 * @brief set the value of attribute 'input_submersion'
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    input_submersion = swString::toLower(text);
    if (!((input_submersion == "resulthydro") ||
          (input_submersion == "resulttime") ||
          (input_submersion == "all")))
    {
        logError() << "  invalid value for parameter 'input_submersion'" << text;
        logError() << "  Allowed values are:";
        logError() << "    'all'\t Compute the submersion maps for both ResultHydro**.txt and ResultTime**.txt files";
        logError() << "    'resulttime'\t  Compute the submersion maps for ResultTime**.txt files only";
        logError() << "    'resulthydro'\t Compute the submersion maps for ResultHydro**.txt files only";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

bool swPostProcessParam::setSubmersionDurations(const std::string& text)
/**
 * @brief set the value of attribute 'submersion_durations'
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    submersion_durations = swString::toLower(text);
    if (!((submersion_durations == "none") ||
          (submersion_durations == "max") ||
          (submersion_durations == "total") ||
          (submersion_durations == "all")))
    {
        logError() << "  invalid value for parameter 'submersion_durations'" << text;
        logError() << "  Allowed values are:";
        logError() << "    'none'\t no computation of submersion durations";
        logError() << "    'max'\t compute the maximum submersion durations";
        logError() << "    'total'\t compute the total submersion durations";
        logError() << "    'all'\t compute both the maximum and total submersion durations";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

bool swPostProcessParam::writeLog()
/**
 * @brief Write all converter parameters into the log folder ()
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    // Open file
    std::ofstream file_parameters(m_io->absolutePath(swFile::log_input_postproc_file));
    if (file_parameters.is_open() == false)
    {
        logError() << "    Can't create file " << swFile::log_input_postproc_file;
        logFatal(__FILE__, __func__, __LINE__);
    }

    // Write parameters
    file_parameters << "compute_max"          << "\t" << std::to_string(compute_max) << std::endl;
    file_parameters << "first_submersion"     << "\t" << std::to_string(first_submersion) << std::endl;
    file_parameters << "h_min"                << "\t" << swString::numberToStr(h_min, 'e', 15) << std::endl;
    file_parameters << "input_format"         << "\t" << input_format << std::endl;
    file_parameters << "input_max"            << "\t" << input_max << std::endl;
    file_parameters << "input_submersion"     << "\t" << input_submersion << std::endl;
    file_parameters << "submersion_durations" << "\t" << submersion_durations << std::endl;

    return EXIT_SUCCESS;
}
