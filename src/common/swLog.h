#pragma once

#include <vector>
#include <fstream>

#include "../common/swMessage.h"

/// define the operator<< for printing message
template <typename T>
inline const swMessage& operator<< (const swMessage& eng, const T& t)
{
    eng.stream() << t << " ";
    *swMessage::outfile << t << " ";
    return eng;
}
template <>
inline const swMessage& operator<< (const swMessage& eng, const std::string& t)
{
    eng.stream() << t.c_str() << " ";
    *swMessage::outfile << t.c_str() << " ";
    return eng;
}
template <>
inline const swMessage& operator<< (const swMessage& eng, const std::vector<std::string> &t)
{
    for (size_t k = 0; k < t.size(); ++k)
    {
        eng.stream() << t[k].c_str() << " ";
        *swMessage::outfile << t[k].c_str() << " ";
    }
    return eng;
}

inline const swMessage& operator<<(const swMessage& eng,std::ostream &(*manip)(std::ostream &))
{
    manip( eng.stream() );
    return eng;
}

template <class T>
inline const swMessage& operator<< (const swMessage& eng, const std::vector<T> &t)
{
    eng.stream() << "[ ";
    *swMessage::outfile << "[ ";
    for (size_t k = 0; k < t.size(); ++k)
    {
        eng.stream() << t[k] << " ";
        *swMessage::outfile << t[k] << " ";
    }
    eng.stream() << "]";
    *swMessage::outfile << "]";

    return eng;
}

template <class T>
inline const swMessage& operator<< (const swMessage& eng, const std::vector<std::vector<T>> &t)
{
    eng.stream() << "[";
    *swMessage::outfile << "[";
    for (size_t k = 0; k < t.size(); ++k)
    {
        eng.stream() << "[";
        *swMessage::outfile << "[";
        for (size_t k_k = 0; k_k < t[k].size(); ++k_k)
        {
            eng.stream() << t[k][k_k] << " ";
            *swMessage::outfile << t[k][k_k] << " ";
        }
        eng.stream() << "]";
        *swMessage::outfile << "]";
    }
    eng.stream() << "]";
    *swMessage::outfile << "]";

    return eng;
}

swMessage logDebug(std::string file, std::string function, int line);
swMessage logEmpt();
swMessage logError();
swMessage logExit();
swMessage logFatal(std::string file, std::string function, int line);
swMessage logInfo();
swMessage logWarn();


inline swMessage logDebug(std::string file, std::string function, int line) {return swMessage(swDebugMsg, file, function, line); } // implementation of logFatal
inline swMessage logEmpt() {return swMessage(swEmptyMsg); } // implementation of logEmty
inline swMessage logError() {return swMessage(swErrorMsg); } // implementation of logError
inline swMessage logExit() {return swMessage(swExitMsg); } // implementation of logExit
inline swMessage logFatal(std::string file, std::string function, int line) {return swMessage(swFatalMsg, file, function, line); } // implementation of logFatal
inline swMessage logInfo() {return swMessage(swInfoMsg); } // implementation of logInfo with new line
inline swMessage logWarn() {return swMessage(swWarningMsg); } // implementation of logWarn with new line
