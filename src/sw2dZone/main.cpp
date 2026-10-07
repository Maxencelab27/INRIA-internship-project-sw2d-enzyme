#include <filesystem>
#include <cstring>

#if defined (WIN32)
#include <windows.h>
#endif

#include "../common/swFile.h"
#include "../common/swIO.h"
#include "../common/swLog.h"
#include "../common/geometry/swMesh.h"
#include "../common/version.h"

#include "swZone.h"
#include "swZonePipeline.h"

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
        std::cerr << "sw2dZone should receive exactly 1 parameter,";
        std::cerr << " received " << argc - 1 << "." << std::endl;
        std::cerr << "  The parameter can be:" << std::endl;
        std::cerr << "    the path of the simulation folder (containing the 'primary_files' and 'input' folders)" << std::endl;
        std::cerr << "    -h or --help : displays the help menu" << std::endl;
        std::cerr << "    -v or --version : show the version of the program" << std::endl;
        exit(EXIT_FAILURE);
    }

    if ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))
    {
        std::cout << std::endl << "DESCRIPTION :" << std::endl;
        std::cout << "  sw2dZone uses the geometry processed by sw2dConverter to generate spatially distributed hydrodynamic parameters over the computational mesh" << std::endl;
        std::cout << std::endl << "PARAMETERS :" << std::endl;
        std::cout << "  sw2dZone should receive exactly 1 parameter" << std::endl;
        std::cout << "  The parameter can be:" << std::endl;
        std::cout << "    the path of the simulation folder (containing the 'primary_files' folder)" << std::endl;
        std::cout << "    -h or --help : show this message" << std::endl;
        std::cout << "    -v or --version : show the version of the program" << std::endl;
        std::cout << std::endl << "CREATED FILES :" << std::endl;
        std::cout << "  A text file in the 'input' folder which name is defined automatically" << std::endl;
        exit(EXIT_SUCCESS);
    }

    if ((strcmp(argv[1], "-v") == 0) || (strcmp(argv[1], "--version") == 0))
    {
        std::cout << "SW2D version " << sw2d_version << std::endl;
        exit(EXIT_SUCCESS);
    }
    //== Descritpion of argv[]

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
    m_io->createEmptyFolders(swFile::log_zone_file);
    swMessage::outfile = new std::ofstream(m_io->absolutePath(swFile::log_zone_file));
    //==== Reset simulation folder
    //== Import from folder name

    //-- Write version and command line to log file
    logInfo() << "Version : " << sw2d_version;
    std::string command_line = argv[0];
    for (int i = 1; i < argc; ++i)
    {
        command_line = command_line + " " + argv[i];
    }
    logInfo() << "Command line:" << command_line;
    //== Write version and command line to log file

    //-- Read the mesh
    logInfo() << "Read the mesh";
    swMesh *m_mesh = new swMesh(m_io);
    logInfo() << "...done";
    if (m_mesh->build() != EXIT_SUCCESS)
    {
        logError() << "Error while building mesh";
        logFatal(__FILE__, __func__, __LINE__);
    }
    //== Read the mesh

    //-- Process
    swZonePipeline *zone_pipeline = new swZonePipeline(m_mesh, m_io);
    if (zone_pipeline->readInput() != EXIT_SUCCESS)
    {
        logError() << "failure during zone_pipeline->readInput()";
        return 1;
    }

    if (zone_pipeline->assignElements() != EXIT_SUCCESS)
    {
        logError() << "failure during zone_pipeline->assignElements()";
        return 1;
    }

    if (zone_pipeline->computeParam() != EXIT_SUCCESS)
    {
        logError() << "failure during zone_pipeline->computeParam()";
        return 1;
    }
    //== Process

    logInfo() << "Done";

    logEmpt() << "Zone ran with " << swMessage::nb_warn << " warning(s)";

    swMessage::finalizeLog();

    return 0;
}

//
// main.cpp ends here
