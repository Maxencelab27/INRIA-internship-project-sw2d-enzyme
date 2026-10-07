#include "swMessage.h"
#include <fstream>
#include "swString.h"

std::ofstream* swMessage::outfile = nullptr;
int swMessage::nb_warn = 0;


swMessage::swMessage(swMessageType type, std::string file, std::string function, int line)
/**
 * @brief constructor of swMessage
 * @param[in] type define the debuging level
 */
{
    if (type == swFatalMsg)
    {
        m_is_fatal = true; // will quit after writing the message
    }

    std::string txt_scr;
    std::string txt_log;

    if (type == swDebugMsg)
    {
        std::string sw2d_file = "src" + swString::splitStr(file, "src")[1];
        txt_scr = "\033[38;5;208m\033[1m[Debug]\033[25m\033[0m " + sw2d_file + " " + function + " " + std::to_string(line) + " ";
        txt_log = "[Debug] " + sw2d_file + " " + function + " " + std::to_string(line) + " ";
    }
    else if (type == swInfoMsg)
    {
        txt_scr = "[Info] ";
        txt_log = "[Info] ";
    }
    else if (type == swWarningMsg)
    {
        txt_scr = "\033[38;5;208m\033[1m[Warn]\033[25m\033[0m ";
        txt_log = "[Warn] ";
        nb_warn++;
    }
    else if (type == swErrorMsg)
    {
        txt_scr = "\033[38;5;196m\033[1m[Err.]\033[25m\033[0m ";
        txt_log = "[Err.] ";
    }
    else if (type == swFatalMsg)
    {
        std::string sw2d_file = "src" + swString::splitStr(file, "src")[1];
        txt_scr = "\033[38;5;13m\033[1;4m[Fat.]\033[25m\033[0m " + sw2d_file + " " + function + " " + std::to_string(line);
        txt_log = "[Fat.] " + sw2d_file + " " + function + " " + std::to_string(line);
    }
    else if (type == swEmptyMsg)
    {
        txt_scr = "";
        txt_log = "       ";
    }
    else if (type == swExitMsg)
    {
        txt_scr = "\033[38;5;10m\033[1m[EXIT]\033[25m\033[0m";
        txt_log = "[EXIT] ";
    }
    stream() << txt_scr;
    *swMessage::outfile << txt_log;
}


swMessage::~swMessage(void) noexcept(false)
/**
 * @brief Destroy the sw Message::sw Message object
 * 
 */
{
    stream() << std::endl;
    *swMessage::outfile << std::endl;
    if (m_is_fatal)
    {
        try
        {
            throw std::runtime_error("Fatal error occured. Quitting.");
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
            exit(1);
        }
    }
}


void swMessage::finalizeLog()
/**
 * @brief close a log message
 * 
 */
{
    outfile->close();
}


std::ostream& swMessage::stream(void) const
/**
 * @brief swMessage::stream
 * @return return the stream in which the message have to be printed
 */
{
    return m_stream;
}
