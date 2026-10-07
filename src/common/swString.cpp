#include "swString.h"

#include <iomanip>
#include <regex>
#include <sstream>

#include "../common/swLog.h"


bool swString::contains(const std::vector<std::string> &s_list, const std::string &s)
/**
 * @brief 
 * 
 * @param[in] s_list 
 * @param[in] s 
 * @return true 
 * @return false 
 */
{
    bool contained = false;
    for (std::string str : s_list)
    {
        if (str == s)
        {
            contained = true;
            break;
        }
    }
    return contained;
}


std::string swString::dirSeparator(void)
/**
 * @brief returns the directory separator (depending on the OS)
 * @return std::string 
 */
{
    #if defined (__linux__)
        return "/";
    #endif
    #if defined (__APPLE__)
        return "/";
    #endif
    #if defined (_WIN32)
        return "\\";
    #endif

    logError() << "Unknown os";
    logFatal(__FILE__, __func__, __LINE__);
}


int swString::indexOf(const std::vector<std::string> list, const std::string str)
/**
 * @brief returns the index of a given string in a string vector
 * @param[in] list vector of strings in which the string 'str' is searched
 * @param[in] str searched string
 * @return the index (ind) of str in the vector (returns -1 is element does not belong to the vector)
 */
{
    size_t i = 0;
    for (std::string string : list)
    {
        if (string == str)
        {
            break;
        }
        ++i;
    }

    if (i > list.size())
    {
        return -1;
    }
    else
    {
        return i;
    }
}


std::string swString::ltrim(const std::string &s) // adapted from https://stackoverflow.com/questions/216823/how-to-trim-an-stdstring
/**
 * @brief removes white space at the beginning of the string
 * @param[in] s string to be cleaned of spaces
 * @return clean string
 */
{
    std::regex regexp(regexWhiteSpace());

    std::string text;
    bool started = false;
    for (size_t i = 0; i < s.size(); ++i)
    {
        std::string substrtxt = s.substr(i, 1);
        if (std::regex_match(substrtxt, regexp) == false)
        {
            started = true;
        }

        if (started == true)
        {
            text += substrtxt;
        }
    }
    return text;
}

template <typename T>

std::string swString::numberToStr(const T& val, char format, int precision)
/**
 * @brief convert a number to string based on the specified format
 */
{
    std::stringstream result;
    switch(format) {
        case 'e': {
            result << std::scientific << std::nouppercase << std::setprecision(precision) << val;
            return result.str();
        }
        case 'f': {
            result << std::fixed << std::setprecision(precision) << val;
            return result.str();
        }
        case 'g': {
            result << std::nouppercase << std::setprecision(precision) << val;
            return result.str();
        }
        default: {
            logError() << "format unknown: " << format;
            return "Format unknown";
        }
    }
}


std::string swString::numberToStr(const double &val, char format, int precision)
/**
 * @brief 
 * @param[in] val 
 * @param[in] format 
 * @param[in] precision 
 * @return std::string 
 */
{
    std::stringstream result;
    switch(format)
    {
        case 'e':
        {
            result << std::scientific << std::nouppercase << std::setprecision(precision) << val;
            return result.str();
        }
        case 'f':
        {
            result << std::fixed << std::setprecision(precision) << val;
            return result.str();
        }
        case 'g':
        {
            result << std::nouppercase << std::setprecision(precision) << val;
            return result.str();
        }
        default:
        {
            logError() << "format unknown: " << format;
            logFatal(__FILE__, __func__, __LINE__);
            return "Format unknown";
        }
    }
}


std::string swString::regexWhiteSpace(void)
/**
 * @brief 
 * @return std::string 
 */
{
    return "\\s+";
}


std::string swString::rtrim(const std::string &s)
/**
 * @brief removes white space at the end of the string
 */
{
    std::string reversed = "";
    for (size_t i = 0; i < s.size(); ++i)
    {
        reversed += s[s.size() - 1 - i];
    }
    std::string reversed_trim = swString::ltrim(reversed);

    std::string text = "";
    for (size_t i = 0; i < reversed_trim.size(); ++i)
    {
        text += reversed_trim[reversed_trim.size() - 1 - i];
    }
    return text;
}


std::vector<std::string> swString::splitStr(const std::string &s, const std::string& pattern)
/**
 * @brief split a string based on a regex separator
 * @param[in] s: string to be splitted
 * @param[in] pattern: regex of the separator
 */{
    std::regex regexp(pattern);
    std::sregex_token_iterator iter(s.begin(), s.end(), regexp, -1);
    std::sregex_token_iterator end;
    return {iter, end};
}


std::vector<std::string> swString::splitStrOnce(const std::string &s, const std::string& pattern, const size_t idx)
/**
 * @brief split a string only once, at the idx separator
 * @param[in] s: string to be splitted
 * @param[in] pattern: regex of the separator
 * @param[in] idx
 */
{
    std::vector<std::string> string_list_temp = splitStr(s, pattern);
    std::vector<std::string> string_list;
    std::string text;
    for (size_t i = 0; i < idx; ++i)
    {
        text += string_list_temp[i];
    }
    string_list.push_back(text);
    text = "";
    for (size_t i = idx; i < string_list_temp.size(); ++i)
    {
        text += string_list_temp[i];
    }
    string_list.push_back(text);

    return string_list;
}


double swString::stringToDouble(const std::string str, bool &conversion_success)
/**
 * @brief convert string to double
 */
{
    double val = 0.0;
    try
    {
        val = std::stod(str);
    }
    catch (std::invalid_argument const& ex)
    {
        logError() << "std::invalid_argument::what(): " << ex.what();
        logError() << __FILE__ << __func__ << __LINE__;
        conversion_success = false;
    }
    return val;
}


long double swString::stringToLongDouble(const std::string str, bool &conversion_success)
/**
 * @brief convert string to long double
 */
{
    long double val = 0.0;
    try
    {
        val = std::stold(str);
    }
    catch (std::invalid_argument const& ex)
    {
        logError() << "std::invalid_argument::what(): " << ex.what();
        logError() << __FILE__ << __func__ << __LINE__;
        conversion_success = false;
    }
    return val;
}


int swString::stringToInt(const std::string str, bool &conversion_success)
/**
 * @brief convert string to int
 */
{
    int val = 0;
    try
    {
        val = std::stoi(str);
    }
    catch (std::invalid_argument const& ex)
    {
        logError() << "std::invalid_argument::what(): " << ex.what();
        logError() << __FILE__ << __func__ << __LINE__;
        conversion_success = false;
    }
    return val;
}


size_t swString::stringToUInt(const std::string str, bool &conversion_success)
/**
 * @brief convert string to size_t
 */
{
    size_t val = 0;
    try
    {
        val = std::stoul(str);
    }
    catch (std::invalid_argument const& ex)
    {
        logError() << "std::invalid_argument::what(): " << ex.what();
        logError() << __FILE__ << __func__ << __LINE__;
        conversion_success = false;
    }
    return val;
}


std::string swString::tab(void)
/**
 * @brief 
 * @return std::string 
 */
{
    return "\t";
}


std::string swString::toLower(const std::string &s)
/**
 * @brief transforms string into lower case
 */
{
    std::string s_to_lower = s;
    std::transform(s_to_lower.begin(), s_to_lower.end(), s_to_lower.begin(), ::tolower);
    return s_to_lower;
}


std::string swString::trim(const std::string &s) // from https://stackoverflow.com/questions/216823/how-to-trim-an-stdstring
/**
 * @brief suppress white space at both ends of the string
 */{
    return swString::ltrim(swString::rtrim(s));
}
