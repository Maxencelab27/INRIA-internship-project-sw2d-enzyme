#include <cstring>
#include <filesystem>
#include <iostream>

#if defined (WIN32)
#include <windows.h>
#endif

#include "../common/swFile.h"
#include "../common/swIO.h"
#include "../common/swLog.h"
#include "../common/version.h"
#include "swPostProcessParam.h"
#include "swPostProcessPipeline.h"

int main(int argc, char **argv)
{
    // -- Handle popups under Windows
#if defined (WIN32)
    SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX | SEM_NOOPENFILEERRORBOX);
#endif
    // == Handle popups under Windows

    //-- Descritpion of argv[]
    if (argc != 2)
    {
        std::cerr << "sw2dPostProcess should receive exactly 1 parameter,";
        std::cerr << " received " << argc << "." << std::endl;
        std::cerr << "  The parameter can be:" << std::endl;
        std::cerr << "    the path of the simulation folder (containing the 'primary_files' folder)" << std::endl;
        std::cerr << "    -h or --help to show this message" << std::endl;
        std::cerr << "    -v or --version to show the version of the program" << std::endl;
        exit(EXIT_FAILURE);
    }

    if ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))
    {
        std::cout << std::endl << "DESCRIPTION :" << std::endl;
        std::cout << "  sw2dPostProcess processes the simulation results to compute things such as maximum values or submersion durations" << std::endl;
        std::cout << std::endl << "PARAMETERS :" << std::endl;
        std::cout << "  sw2dPostProcess should receive exactly 1 parameter" << std::endl;
        std::cout << "  The parameter can be:" << std::endl;
        std::cout << "    the path of the simulation folder (containing the 'input' and 'output' folders)" << std::endl;
        std::cout << "    -h or --help to show this message" << std::endl;
        std::cout << "    -v or --version to show the version of the program" << std::endl;
        std::cout << std::endl << "CREATED FILES :" << std::endl;
        std::cout << "  ['output/ResultHydro_Max.txt' and 'output/ResultHydro_DateOfMax.txt'] and/or SubmersionDurations.txt" << std::endl;
        exit(EXIT_SUCCESS);
    }
    //== Descritpion of argv[]

    if ((strcmp(argv[1], "-v") == 0) || (strcmp(argv[1], "--version") == 0))
    {
        std::cout << "SW2D version " << sw2d_version << std::endl;
        exit(EXIT_SUCCESS);
    }

    //-- Import from folder name
    const std::string user_defined_folder = argv[1];
    //---- Check existing folder
    if (std::filesystem::exists(user_defined_folder) == false)
    {
        std::cerr << "Folder " << user_defined_folder  << " does not exist" << std::endl;
        std::cerr << "Use option -h or --help to display option menu" << std::endl;
        return EXIT_FAILURE;
    }
    //==== Check existing folder

    //------ Reset simulation folder
    swIO *m_io = new swIO(user_defined_folder);
    m_io->createEmptyFolders(swFile::log_postprocess_file);
    swMessage::outfile = new std::ofstream(m_io->absolutePath(swFile::log_postprocess_file));
    //==== Reset simulation folder

    //== Import from folder name

    //---- Write version and command line to log file
    logInfo() << "Version : " << sw2d_version;
    std::string command_line = argv[0];
    for (int i = 1; i < argc; ++i)
    {
        command_line = command_line + " " + argv[i];
    }
    logInfo() << "Command line:" << command_line;
    //==== Write version and command line to log file

    //Initialize the pointer to swPostProcessParam class
    swPostProcessParam *m_param = new swPostProcessParam(m_io);

    //-- reading parameters
    logInfo() << "Importing post process parameters from file:" << swFile::input_postprocess_file;
    if (m_param->readInputFile() != EXIT_SUCCESS)
    {
        logError() << "Failure while reading port process parameters from file" << swFile::input_postprocess_file;
        logFatal(__FILE__, __func__, __LINE__);
    }
    logInfo() << "...done";
    //== reading parameters

    //-- check parameters consistency
    //---- check consistency
    /*logInfo() << "Check input parameters consistency";
    if (m_param->checkConsistency() == false)
    {
        logError() << "  incoherence detected.";
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }*/
    //==== check consistency
    //== check parameters consistency

    //-- write the input parameters
    logInfo() << "Writting post process parameters into file:" << m_io->absolutePath(swFile::log_input_postproc_file);
    if (m_param->writeLog() != EXIT_SUCCESS)
    {
        logError() << "Failure while writting postprocess parameters into file" << swFile::log_input_postproc_file;
        logFatal(__FILE__, __func__, __LINE__);
    }
    logInfo() << "...done";
    //== write the input parameters

    swPostProcessPipeline *m_postproc = new swPostProcessPipeline(m_io, m_param);

    //-- Compute max
    if (m_param->getComputeMax() == true)
    {
        if(m_postproc->computeMax() != EXIT_SUCCESS)
        {
            logError() << "Failure while computing max";
            logFatal(__FILE__, __func__, __LINE__);
        }
    }
    //== Compute max

    //-- Compute submersion duration
    if ((m_param->getFirstSubmersion()) || (m_param->getSubmersionDurations() != "none"))
    {
        if (m_postproc->computeSubmersions() != EXIT_SUCCESS)
        {
            logError() << "Failure while computing submersion durations";
            logFatal(__FILE__, __func__, __LINE__);
        }
    }
    //== Compute submersion duration

	return 0;
}

// main.cpp ends here
