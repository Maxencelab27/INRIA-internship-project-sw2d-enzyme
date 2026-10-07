#pragma once

#include <string>

class swIO;

class swPorosityMapperParam
{
public :
    swPorosityMapperParam(swIO *io);
    ~swPorosityMapperParam(void) = default;

    // General methods
    bool checkConsistency();
    bool readInputFile(std::string filepath);
    bool writeLog(std::string filepath);

    // Getters
    std::string getInput();
    std::string getModel();
private:
    //-- Local attributes
    swIO *m_io = nullptr;              /**< pointer to the swIO object to handle files IO */
    //== Local attributes

    //-- User-defined attributes
    std::string input = "footprint";   /**< format for the buildings input */
    std::string model = "none";        /**< name of the porosity model */
    //== User-defined attributes

};
