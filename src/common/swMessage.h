#pragma once

//#include <swCoreExport>
#include <iostream> // for std::cerr
// #include <regex>
// #include <vector>
// #include <array>
//#include <string>
//#include <fstream>

/**
 * @class swMessage
 * @brief Dedicated class to provide various levels of messages (info, warning, error, fatal)
 */
enum swMessageType {swDebugMsg, swEmptyMsg, swExitMsg, swInfoMsg, swWarningMsg, swErrorMsg, swFatalMsg };

class swMessage
{
public:
    //-- General methods
    swMessage(swMessageType type, std::string file = "", std::string function = "", int line = -1);
    virtual ~swMessage(void) noexcept(false);

    static void finalizeLog();
    std::ostream& stream(void) const;
    //== General methods

    //-- Attributes
    static int nb_warn;
    static std::ofstream* outfile;
    //== Attributes

private:
    //-- Attributes
    bool m_is_fatal = false; // exit on fatal error messages
    std::ostream& m_stream = std::cerr; // instanciate the output stream.
    //== Attributes
};
