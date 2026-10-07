#include <cstring>
#include <filesystem>
#include <iostream>

#if defined (WIN32)
#include <windows.h>
#endif

#include "../common/swFile.h"
#include "../common/swIO.h"
#include "../common/swLog.h"
#include "../common/geometry/swMesh.h"
#include "../common/version.h"

#include "swPorosityMapperParam.h"
#include "swPorosityMapperPipeline.h"

int main(int argc, char **argv)
{
    // -- Handle popups under Windows
#if defined (WIN32)
    SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX | SEM_NOOPENFILEERRORBOX);
#endif
    // == Handle popups under Windows

    //-- Descritpion of argv[]
    if (!(argc == 2))
    {
        std::cerr << "sw2dConverter should receive 1 parameter";
        std::cerr << " and received " << argc - 1 << std::endl;
        std::cerr << "  The parameter can be:" << std::endl;
        std::cerr << "    the path of the simulation folder (containing the 'primary_files' folder)" << std::endl;
        std::cerr << "    -h or --help : displays the help menu" << std::endl;
        std::cerr << "    -v or --version : show the version of the program" << std::endl;
        exit(EXIT_FAILURE);
    }

    if ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))
    {
        std::cout << std::endl << "DESCRIPTION :" << std::endl;
        std::cout << "  sw2dConverter processes the mesh so as to derive all the necessary numerical parameters prior to the computations" << std::endl;
        std::cout << std::endl << "PARAMETER :" << std::endl;
        std::cout << "  sw2dConverter should receive exactly 1 parameter" << std::endl;
        std::cout << "  The parameter can be::" << std::endl;
        std::cout << "    the path of the simulation folder (containing the 'primary_files' folder)" << std::endl;
        std::cout << "    -h or --help : show this message" << std::endl;
        std::cout << "    -v or --version : show the version of the program" << std::endl;
        std::cout << std::endl << "CREATED FILES :" << std::endl;
        std::cout << "  " << swFile::input_geo_file << std::endl;
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
    m_io->createEmptyFolders(swFile::log_porositymapper_file);
    swMessage::outfile = new std::ofstream(m_io->absolutePath(swFile::log_porositymapper_file));
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

    //Initialize the pointer to swPorosityMapperParam object
    swPorosityMapperParam *m_param = new swPorosityMapperParam(m_io);

    //-- reading parameters
    logInfo() << "Importing porosity mapper parameters from file:" << swFile::input_porosity_mapper_file;
    if (m_param->readInputFile(swFile::input_porosity_mapper_file) != EXIT_SUCCESS)
    {
        logError() << "Failure during swConverterParam::read";
        logFatal(__FILE__, __func__, __LINE__);
    }
    logInfo() << "...done";
    //== reading parameters

    //-- check parameters consistency
    //---- check consistency
    logInfo() << "Check input parameters consistency";
    if (m_param->checkConsistency() == false)
    {
        logError() << "  incoherence detected.";
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    //==== check consistency
    //== check parameters consistency

    //-- write the input parameters
    logInfo() << "Writting porosity mapper parameters into file:" << m_io->absolutePath(swFile::log_input_porosity_mapper_file);
    if (m_param->writeLog(swFile::log_input_porosity_mapper_file) != EXIT_SUCCESS)
    {
        logError() << "Failure during swConverterParam::write";
        logFatal(__FILE__, __func__, __LINE__);
    }
    logInfo() << "...done";
    //== write the input parameters

    swMesh *m_mesh = new swMesh(m_io);

    //-- import the mesh
    logInfo() << "Import mesh";
    if (m_mesh->build() != EXIT_SUCCESS)
    {
        logError() << "Error while building mesh";
        logFatal(__FILE__, __func__, __LINE__);
    }
    logInfo() << "...done";
    //== import the mesh

    //-- initialise the pipeline
    swPorosityMapperPipeline *porositymapper_pipeline = new swPorosityMapperPipeline(m_io, m_param, m_mesh);
    porositymapper_pipeline->init();
    //== initialise the pipeline

    //-- map porosity
    if (porositymapper_pipeline->mapPorosity(m_param->getModel()) == EXIT_FAILURE)
    {
        logError() << "failure during porositymapper_pipeline->computePorosity()";
        return 1;
    }
    //== map porosity

    logEmpt() << "Porosity mapper ran with " << swMessage::nb_warn << " warning(s)";
    swMessage::finalizeLog();

    return EXIT_SUCCESS;

	return 0;
}

// main.cpp ends here
