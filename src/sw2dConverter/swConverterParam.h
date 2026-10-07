#include <string>

class swIO;

class swConverterParam
{
public :
    swConverterParam(swIO *io);
    ~swConverterParam(void) = default;

    bool checkConsistency();
    bool readInputFile();
    bool writeLog();

    bool correctTopography();
    double getAMin();
    std::string getBcFileName();
    int getDefaultBc();
    std::string getMeshFileName();
    bool getOverWrite();
    std::string getTopoFileName();
    std::string getTopoOption();
    bool getWriteLogs();

private:
    //-- Local attributes
    bool default_bc_defined = false;   /**< true if 'default_bc' has been user defined */
    swIO *m_io = nullptr;              /**< pointer to the swIO object to handle files IO */
    bool topo_correction = false;      /**< flag to activate topo_correction */
    //== Local attributes

    //-- User-defined attributes
    double a_min = -1;                  /**< threshold for 'small cells' */
    std::string bc_file = "";           /**< bc file name in primary_files directory */
    int default_bc = 0;                 /**< default value for edges bc number */
    std::string mesh_file = "";         /**< mesh file name in primary_files directory */
    bool overwrite = false;             /**< flag to overwrite previous geometry file */
    std::string topo_file = "";         /**< topo file name in primary_files directory */
    std::string topo_option = "";       /**< option for topography correction */
    bool write_logs = true;             /**< flag to write mesh log files */
    //== User-defined attributes

};
