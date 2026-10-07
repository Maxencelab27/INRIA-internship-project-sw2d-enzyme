#pragma once
#include <string>
#include <vector>

class swIO
{
public:
    ///-- General Methods
    swIO(const std::string &user_defined_folder);
    ~swIO(void);
    std::string absolutePath(const std::string &file);
    void computeFileNameFromTime(double const& t_time, std::string& file_name);
    void computeFileNameFromTime(double const& t_time, std::string label, std::string& file_name);
    void createEmptyFolders(const std::string &log_file);
    void createEmptyLogFolder(const std::string &log_file);
    void createEmptyOutputFolder(bool remove_output);
    static bool createTxtFileMapDecimal(std::string const& file_name, std::vector<std::string> const& var_names,
                                        std::vector<std::vector<double>> const& var_values, size_t precision);
    
    static bool createFileMapString(std::string const& file_name, std::vector<std::string> const& var_names,
                                    std::vector<std::vector<std::string>> const& var_values);
    static bool fillBinFileMapDecimal(const std::string &filename, double time, std::vector<std::vector<double>> const& var_values);
    static bool initBinFileMapDecimal(const std::string &filename, const std::string &geo_filename, const std::vector<std::string> &var_names,
                                      const double t_0, const double t_max, const double dt_map);
    bool readBoundaryRatingCurveFile(const std::string& filename, std::vector<std::vector<std::string>>& rating_curve_data);
    bool readBoundaryTimeSeriesFile(const std::string& filename,
                                    std::vector<int> &bc_timeseries_codes,
                                    std::vector<std::string> &bc_timeseries_type,
                                    std::vector<double> &timebc,
                                    std::vector<std::vector<double>> &val_bc);
    static bool readTimeSeriesScalarFile(const std::string& filename,
                                         std::vector<double>& time_vector,
                                         std::vector<std::vector<double>>& val_table);
    static bool readTimeSeries2dVectorFile(const std::string& filename,
                                         std::vector<double>& time_vector,
                                         std::vector<std::vector<double>>& val_1_table,
                                         std::vector<std::vector<double>>& val_2_table);
    static void readVectorParameterFile(std::string const& file_name, size_t const& n_params,
                                        size_t const& n_records, std::vector<std::string>& param_name,
                                        std::vector<std::vector<double>>& param_array); 
    static double timeFromFileName(std::string filename);
    static void errorReadingLine(std::ifstream& file,
                                    const std::string& filename,
                                    const std::string& line_txt,
                                    int line_input_file, int line_code,
                                    const std::string& function_name,
                                    const std::string &cpp_name);
    static std::string readLine(std::ifstream& file);
    static std::vector<std::string> stringSplit(std::string const& str);
    static bool writeTxtLine(std::ofstream& output_stream, std::string file_line);
    //== General methods

public:
    ///-- Attributes
    std::string project_folder;      /**< simulation directory path (relative to where the executable was called), terminated by a folder separator **/
    //== Attributes
};
