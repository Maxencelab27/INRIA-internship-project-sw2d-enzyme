#pragma once

#include <string>
#include <vector>

namespace swString
{
    ///-- General methods
    // separator chars
    std::string regexWhiteSpace(void);
    std::string tab(void);
    std::string dirSeparator(void);
    
    /// string list
    bool contains(const std::vector<std::string> &s_list, const std::string &s);
    int indexOf(const std::vector<std::string> list, const std::string str);
    std::vector<std::string> splitStr(const std::string &s,const std::string &pattern);
    std::vector<std::string> splitStrOnce(const std::string &s,const std::string &pattern, const size_t idx);

    /// conversion
    template <typename T>
    std::string numberToStr(const T& val, char format, int precision);
    std::string numberToStr(const double& val, char format, int precision);
    double stringToDouble(const std::string str, bool &conversion_success);
    long double stringToLongDouble(const std::string str, bool &conversion_success);
    int stringToInt(const std::string str, bool &conversion_success);
    size_t stringToUInt(const std::string str, bool &conversion_success);

    /// remove white space at the ends of string
    std::string ltrim(const std::string &s);
    std::string rtrim(const std::string &s);
    std::string trim(const std::string &s);

    /// get string to lower
    std::string toLower(const std::string &s);
    //== General methods
}

inline const std::vector<std::string>& operator+ (std::vector<std::string>& list1, const std::vector<std::string>& list2)
{
    std::vector<std::string> list3 = list1;
    for (std::string text:list2)
    {
        list1.push_back(text);
    }
    return list1;
}
