#pragma once

#include <string>
#include <vector>

#include "swPostProcessParam.h"
#include "../common/swIO.h"

class swPostProcessPipeline
{
public:
    //-- General methods
    swPostProcessPipeline(swIO *io, swPostProcessParam *param);
    ~swPostProcessPipeline() = default;

    bool computeMax();
    bool computeSubmersions();
    //== General methods

private:
    //-- General methods
    bool importResults(std::string format_option, std::string files_option);
    bool importBinResults();
    bool importTxtResults(std::string files_option);
    //== General methods

    //--Attributes
    swPostProcessParam *m_param;
    swIO *m_io;

    bool result_hydro_imported = false;    /** < indicating whether ResultHydro maps have been imported or not */
    bool result_time_imported = false;     /** < indicating whether ResultTime maps have been imported or not */
    bool result_bin_imported = false;      /** < indicating whether binary results have been imported or not */

    std::vector<double> x_cells;                            /** < x-coordinates of cells in results */
    std::vector<double> y_cells;                            /** < y-coordinates of cells in results */
    std::vector<double> zb_cells;                           /** < bottom elevation of cells in results */
    std::vector<std::string> fields;                        /** < list of results fields */
    std::vector<double> times;                              /** < list of results times */
    std::vector<std::vector<std::vector<double>>> values;   /** < list of results nb_times x nb_fields x nb_cells */
    //==Attributes
};
