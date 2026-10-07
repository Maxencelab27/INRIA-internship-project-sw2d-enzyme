#include "swSolverParam.h"

#include <filesystem>

#include "../../common/swFile.h"
#include "../../common/swIO.h"
#include "../../common/swLog.h"


swSolverParam::swSolverParam(swIO *io)
/**
 * @brief create solverParam object
 * 
 */
{
    m_io = io;
    logInfo() << "Simulation folder: " << m_io->project_folder << std::endl;
}


bool swSolverParam::checkConsistency()
/**
 * @brief checks consistency conditions between input parameters
 * @return bool of exit_success/failure
 */
{
    if (t_max < t_0)
    {
        logError() << "   Maximum simulation time must be larger than initial time: t_max > t_0";
        logFatal(__FILE__, __func__, __LINE__);
    }
    if ((divergence_correction == "limit_flux") && (n_it_div_max <= 0))
    {
        logError() << "Max. number of iteration expected (divergence correction is activated)";
        logError() << "   Please correct parameter 'n_it_div_max'";
        logFatal(__FILE__, __func__, __LINE__);
    }
    if (dt_map > t_max - t_0)
    {
        logWarn() << "   Time step for map storage greater than Tmax => no maps will be stored";
    }
    if (dt_time_series > t_max - t_0)
    {
        logWarn() << "   Time step for probes storage greater than Tmax => no probes will be stored";
    }
    if (cfl_max > 1.)
    {
        /** @todo is it 1 for any model? If it is model dependant,
         *  we need to make the verification using overloaded function of swHyp **/
        logWarn() << "   Courant number defined value: CFLmax=" << cfl_max;
        logWarn() << "     CFLmax has to be smaller than 1 for stability constraint";
        logWarn() << "     CFLmax has been reduced to 1";
        cfl_max = 1;
    }
    return EXIT_SUCCESS;
}


double swSolverParam::getCFLMax()
/**
 * @brief cfl_max getter function
 * @return double cfl_max 
 */
{
    return cfl_max;
}


std::vector<std::string> swSolverParam::getBiblioReferenceList()
/**
 * @brief Bibliography references getter function
 * @return std::vector<std::string> biblio_reference_list
 */
{
    return biblio_reference_list;
}


double swSolverParam::getDdpChiPhiMax()
/**
 * @brief phi ratio maximum value in DDP model
 * @return double chi_phi_max 
 */
{
    return ddp_chi_phi_max;
}


double swSolverParam::getDdpChiThetaMax()
/**
 * @brief theta ratio maximum value in DDP model
 * @return double chi_theta_max 
 */
{
    return ddp_chi_theta_max;
}


std::string swSolverParam::getDivergenceCorrection()
/**
 * @brief Divergence correction getter function
 * @return std::string divergence_correction 
 */
{
    return divergence_correction;
}


double swSolverParam::getDtMap()
/**
 * @brief dt_map getter function
 * @return double dt_map 
 */
{
    return dt_map;
}


double swSolverParam::getDtMax()
/**
 * @brief dt_max getter function
 * @return double dt_max 
 */
{
    return dt_max;
}


double swSolverParam::getDtTimeSeries()
/**
 * @brief dt_time_series getter function
 * @return double dt_time_series 
 */
{
    return dt_time_series;
}

std::vector<std::pair<unsigned long,std::string>>& swSolverParam::getExtractionTimes(void)
/**
 * @brief extraction_times getter function
 * @return std::vector<std::pair<unsigned long,std::string>>& extraction_times 
 */
{
    return extraction_times;
}


double swSolverParam::getFrMax()
/**
 * @brief Fr_max getter function
 * @return double Fr_max 
 */
{
    return Fr_max;
}


double swSolverParam::getGravity()
/**
 * @brief Gravity getter function
 * @return double g=gravity 
 */
{
    return g;
}

bool swSolverParam::setHasQsection(const std::string &text)
/**
 * @brief Sets the boolean q_section (default false)
 * @param[in] text containing the value
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    int val = swString::stringToInt(text, conversion_success);
    if ((!conversion_success) || !((val == 0) || (val == 1)))
    {
        logInfo() << "  unallowed value for parameter 'q_section'" << text;
        logInfo() << "    allowed values for parameter 'q_section'";
        logInfo() << "    - 0\tno q sections";
        logInfo() << "    - 1\tq sections are computed";
        return EXIT_FAILURE;
    }
    if (val == 0)
    {
        has_qsection = false;
    }
    if (val == 1)
    {
        has_qsection = true;
    }
    return EXIT_SUCCESS;
}


double swSolverParam::getHMin()
/**
 * @brief h_min getter function
 * @return double h_min 
 */
{
    return h_min;
}


double swSolverParam::getHMinDdp()
/**
 * @brief h_min_ddp getter function
 * @return double h_min_ddp
 */
{
    return h_min_ddp;
}


bool swSolverParam::getHasQsection()
/**
 * @brief has_qsection getter function
 * @return bool has_qsection 
 */
{
    return has_qsection;
}


bool swSolverParam::getHotstart()
/**
 * @brief hotstart getter function
 * @return bool hotstart 
 */
{
    return hotstart;
}


bool swSolverParam::getLogWrite()
/**
 * @brief log_write getter function
 * @return bool log_write 
 */
{
    return log_write;
}


std::string swSolverParam::getModelBuildingsExchange()
/**
 * @brief model_buildings_exchange getter function
 * @return std::string model_buildings_exchange 
 */
{
    return model_buildings_exchange;
}


std::string swSolverParam::getModelFriction()
/**
 * @brief model_friction getter function
 * @return std::string model_friction 
 */
{
    return model_friction;
}


std::string swSolverParam::getModelHyp()
/**
 * @brief model_hyp getter function
 * @return std::string model_hyp 
 */
{
    return model_hyp;
}


std::string swSolverParam::getModelInfiltration()
/**
 * @brief model_infiltration getter function
 * @return std::string model_infiltration 
 */
{
    return model_infiltration;
}


std::string swSolverParam::getModelLogJam()
/**
 * @brief model_logjam getter function
 * @return std::string model_logjam 
 */
{
    return model_logjam;
}


std::string swSolverParam::getModelMomentumDiffusion()
/**
 * @brief model_momentum_diffusion getter function
 * @return std::string model_momentum_diffusion 
 */
{
    return model_momentum_diffusion;
}


std::string swSolverParam::getModelMomentumDissipation()
/**
 * @brief momentum_dissipation getter function
 * @return bool momentum_dissipation 
 */
{
    return model_momentum_dissipation;
}


std::string swSolverParam::getModelParab()
/**
 * @brief model_parab getter function
 * @return std::string model_parab 
 */
{
    return model_parab;
}


std::string swSolverParam::getModelPorosityRainfall()
/**
 * @brief model_rainfall_porosity getter function
 * @return std::string model_rainfall_porosity 
 */
{
    return model_porosity_rainfall;
}


std::string swSolverParam::getModelRunoff()
/**
 * @brief model_runoff getter function
 * @return std::string model_runoff 
 */
{
    return model_runoff;
}


std::string swSolverParam::getModelSingularHeadLoss()
/**
 * @brief model_singular_head_loss getter function
 * @return std::string model_singular_head_loss 
 */
{
    return model_singular_head_loss;
}


std::string swSolverParam::getModelTransportReaction()
/**
 * @brief model_transport_reaction getter function
 * @return std::string model_transport_reaction  
 */
{
    return model_transport_reaction;
}


std::string swSolverParam::getModelWind()
/**
 * @brief model_wind getter function
 * @return std::string model_wind  
 */
{
    return model_wind;
}


double swSolverParam::getMomentumDiffusionCoefficient()
/**
 * @brief momentum_diffusion_coefficient getter function
 * @return double momentum_diffusion_coefficient 
 */
{
    return momentum_diffusion_coefficient;
}


size_t swSolverParam::getMomentumDiffMaxIter()
/**
 * @brief momentum_diff_max_iter getter function
 * @return size_t momentum_diff_max_iter 
 */
{
    return momentum_diff_max_iter;
}


double swSolverParam::getMomentumDiffConvergenceEps()
/**
 * @brief momentum_diff_convergence_eps getter function
 * @return double momentum_diff_convergence_eps 
 */
{
    return momentum_diff_convergence_eps;
}

size_t swSolverParam::getNDecimal()
/**
 * @brief n_decimal getter function
 * @return size_t n_decimal 
 */
{
    return n_decimal;
}

size_t swSolverParam::getNItDivMax()
/**
 * @brief n_it_div_max getter function
 * @return size_t n_it_div_max 
 */
{
    return n_it_div_max;
}


size_t swSolverParam::getNPorTabMax()
/**
 * @brief n_por_tab_max getter function
 * @return size_t n_por_tab_max 
 */
{
    return n_por_tab_max;
}


std::string swSolverParam::getNegativeDepthDetection()
/**
 * @brief negative_depth_detection getter function
 * @return std::string negative_depth_detection 
 */
{
    return negative_depth_detection;
}


std::string swSolverParam::getNumScheme()
/**
 * @brief num_scheme getter function
 * @return std::string num_scheme
 */
{
    return num_scheme;
}


std::string swSolverParam::getOutputFormat()
/**
 * @brief output_format getter function
 * @return std::string output_format 
 */
{
    return output_format;
}


double swSolverParam::getPhiMin()
/**
 * @brief phi_min getter function
 * @return double phi_min
 */
{
    return phi_min;
}



std::vector<std::string> swSolverParam::getProbeNames()
/**
 * @brief probe_names getter function
 * @return std::vector<std::string> probe_names 
 */
{
    return probe_names;
}


std::vector<double> swSolverParam::getProbeX()
/**
 * @brief probe_x getter function
 * @return std::vector<double> probe_x 
 */
{
    return probe_x;
}


std::vector<double> swSolverParam::getProbeY()
/**
 * @brief probe_y getter function
 * @return std::vector<double> probe_y 
 */
{
    return probe_y;
}


std::string swSolverParam::getStoreFlag()
/**
 * @brief store_flag getter function
 * @return std::string store_flag 
 */
{
    return store_flag;
}


double swSolverParam::getT0()
/**
 * @brief t_0 getter function
 * @return double t_0 
 */
{
    return t_0;
}


double swSolverParam::getTMax()
/**
 * @brief t_max getter function
 * @return double t_max 
 */
{
    return t_max;
}


std::string swSolverParam::getTopoDiscretization()
/**
 * @brief topo_discretization getter function
 * @return std::string  
 */
{
    return topo_discretization;
}


double swSolverParam::getTransportDiffusionCoefficient()
/**
 * @brief transport_diffusion_coefficient getter function
 * @return double transport_diffusion_coefficient 
 */
{
    return transport_diffusion_coefficient;
}


double swSolverParam::getVMax()
/**
 * @brief v_max getter function
 * @return double v_max 
 */
{
    return v_max;
}


bool swSolverParam::readExtractionTime()
/**
 * @brief read additional extraction times (for which a result map is expected to be written)
 * @return bool of exit_success/failure
 */
{
    // empty the extraction_times variable
    extraction_times.clear();

    // set path
    const std::string filename = m_io->absolutePath(swFile::input_extraction_time_file);

    logInfo() << "Check for (not mandatory) Extra Output Times from:" << filename;
    if (std::filesystem::exists(filename) == true)
    {
        // open file
        std::ifstream file(filename);
        if (file.is_open() == false)
        {
            logError() << "Unable to open the input file: " << swFile::input_extraction_time_file;
            logFatal(__FILE__, __func__, __LINE__);
            return EXIT_FAILURE;
        }

        // file parsing
        size_t line_nb = 0;  // line increment is done when starting the loop

        while (file.eof() == false)
        {
            size_t time;
            std::string label;
            bool conversion_success = true;
            std::string line_txt = swIO::readLine(file);
            line_nb++;

            if ((line_txt.empty() == false) && (line_txt.at(0) != '#'))
            {
                std::vector<std::string> list = swString::splitStrOnce(swString::trim(line_txt), swString::regexWhiteSpace(), 1);
                if (list.size() != 2)
                {
                    logError() << "formatting error at line" << line_nb;
                    return EXIT_FAILURE;
                }
                time = swString::stringToUInt(list[0], conversion_success);
                if (conversion_success == false)
                {
                    logError() << "line is misformatted: format conversion";
                    swIO::errorReadingLine(file, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
                }
                label = list[1];

                // insert data into list
                if (extraction_times.size() == 0)
                {
                    extraction_times.push_back(std::pair(time, label));
                }
                else
                {
                    // sort it in place
                    bool greater = true;
                    for (size_t i = 0 ; i < extraction_times.size(); ++i)
                    {
                        size_t t = extraction_times[i].first;

                        if (time == t)
                        {
                            // ignore multiple entry for the same time
                            greater = false;
                            logWarn() << "time" << time << "duplicated at line" << line_nb << "(IGNORED)";
                            break;
                        }

                        if (time < t)
                        {
                            greater = false;
                            extraction_times.insert(extraction_times.begin() + i, std::pair(time, label));
                            break;
                        }
                    }

                    if (greater == true)
                    {
                        extraction_times.push_back(std::pair(time, label));
                        continue;
                    }
                }
            }
        }
        file.close();
        logInfo() << "...done";
    }
    else
    {
        logInfo() << "... no Extra Output Times file provided";
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::readInputFile()
/**
 * @brief Reads the input.solv file and sets the corresponding solver parameters
 * @return bool of exit_success/failure
 * @par Modified attributes
 *  - most of this class attributes
 */
{
    std::ifstream file;
    file.open(m_io->absolutePath(swFile::input_solver_file));
    if (file.is_open() == false)
    {
        logError() << "    Unable to open:" << swFile::input_solver_file;
        logFatal(__FILE__, __func__, __LINE__);
    }
    else
    {
        logInfo() << "Importing parameters from file:"
                  << swFile::input_solver_file;

        size_t line_nb = 1;
        size_t arret = 0;
        while (!file.eof())
        {
            std::string line_txt = swIO::readLine(file);

            if (!(line_txt.empty() || (line_txt[0] == '#') || (line_txt[0]== '=')))
            {
                std::vector<std::string> line_split = swString::splitStr(line_txt, swString::regexWhiteSpace());
                if (line_split.size() != 2)
                {
                    logError() << "line is misformatted";
                    return EXIT_FAILURE;
                }
                else
                {
                    //-- testing the allowed parameters
                    std::string parameter_name = swString::toLower(line_split[0]);
                    if (parameter_name == "cfl_max")
                    {
                        arret += setCFLMax(line_split[1]);
                    }
                    else if (parameter_name == "contact_restitution")
                    {
                        arret += setCoeffRestitution(line_split[1]);
                    }
                    else if (parameter_name == "contact_ext_margin")
                    {
                        arret += setContactExtMargin(line_split[1]);
                    }
                    else if (parameter_name == "ddp_chi_phi_max")
                    {
                        arret += setDdpChiPhiMax(line_split[1]);
                    }
                    else if (parameter_name == "ddp_chi_theta_max")
                    {
                        arret += setDdpChiThetaMax(line_split[1]);
                    }
                    else if (parameter_name == "ddp_ztab0_log")
                    {
                        arret += setDdpZtab0Log(line_split[1]);
                    }
                    else if (parameter_name == "div_corr")
                    {
                        arret += setDivergenceCorrection(line_split[1]);
                    }
                    else if (parameter_name == "dt_map")
                    {
                        arret += setDtMap(line_split[1]);
                    }
                    else if (parameter_name == "dt_max")
                    {
                        arret += setDtMax(line_split[1]);
                    }
                    else if (parameter_name == "dt_time_series")
                    {
                        arret += setDtTimeSeries(line_split[1]);
                    }
                    else if (parameter_name == "fr_max")
                    {
                        arret += setFrMax(line_split[1]);
                    }
                    else if (parameter_name == "gravity")
                    {
                        arret += setGravity(line_split[1]);
                    }
                    else if (parameter_name == "h_min")
                    {
                        arret += setHmin(line_split[1]);
                    }
                    else if (parameter_name == "hotstart_time")
                    {
                        arret += setHotstart(line_split[1]);
                    }
                    else if (parameter_name == "log_write")
                    {
                        arret += setLogWrite(line_split[1]);
                    }
                    else if (parameter_name == "model_building")
                    {
                        arret += setModelBuilding(line_split[1]);
                    }
                    else if (parameter_name == "model_friction")
                    {
                        arret += setModelFriction(line_split[1]);
                    }
                    else if (parameter_name == "model_hyperb")
                    {
                        arret += setModelHyperb(line_split[1]);
                    }
                    else if (parameter_name == "model_infiltration")
                    {
                        arret += setModelInfiltration(line_split[1]);
                    }
                    else if (parameter_name == "model_logjam")
                    {
                        arret += setModelLogJam(line_split[1]);
                    }
                    else if (parameter_name == "model_porosity_rainfall")
                    {
                        arret += setModelPorosityRainfall(line_split[1]);
                    }
                    else if (parameter_name == "model_runoff")
                    {
                        arret += setModelRunoff(line_split[1]);
                    }
                    else if (parameter_name == "model_singular_head_loss")
                    {
                        arret += setModelSingularHeadLoss(line_split[1]);
                    }
                    else if (parameter_name == "model_transport_reaction")
                    {
                        arret += setModelTransportReaction(line_split[1]);
                    }
                    else if (parameter_name == "model_transport_diffusion")
                    {
                        arret += setModelTransportDiffusion(line_split[1]);
                    }
                    else if (parameter_name == "model_wind")
                    {
                        arret += setModelWind(line_split[1]);
                    }
                    else if (parameter_name == "model_momentum_diffusion")
                    {
                        arret += setMomentumDiffusion(line_split[1]);
                    }
                    else if (parameter_name == "model_momentum_dissipation")
                    {
                        arret += setMomentumDissipation(line_split[1]);
                    }
                    else if (parameter_name == "momentum_diffusion_coefficient")
                    {
                        arret += setMomentumDiffusionCoefficient(line_split[1]);
                    }
                    else if (parameter_name == "negative_depth_detection")
                    {
                        arret += setNegativeDepthDetection(line_split[1]);
                    }
                    else if (parameter_name == "n_decimal")
                    {
                        arret += setNDecimal(line_split[1]);
                    }
                    else if (parameter_name == "n_it_div_max")
                    {
                        arret += setNitDivMax(line_split[1]);
                    }
                    else if (parameter_name == "n_por_tab_max")
                    {
                        arret += setNPorTabMax(line_split[1]);
                    }
                    else if (parameter_name == "num_scheme")
                    {
                        arret += setNumScheme(line_split[1]);
                    }
                    else if (parameter_name == "output_format")
                    {
                        arret += setOutputFormat(line_split[1]);
                    }
                    else if (parameter_name == "phi_min")
                    {
                        arret += setPhiMin(line_split[1]);
                    }
                    else if (parameter_name == "q_section")
                    {
                        arret += setHasQsection(line_split[1]);
                    }
                    else if (parameter_name == "read_wall_from_file")
                    {
                        arret += setReadWallFromFile(line_split[1]);
                    }
                    else if (parameter_name == "store_flag")
                    {
                        arret += setStoreFlag(line_split[1]);
                    }
                    else if (parameter_name == "t_0")
                    {
                        arret += setT0(line_split[1]);
                    }
                    else if (parameter_name == "t_max")
                    {
                        arret += setTmax(line_split[1]);
                    }
                    else if (parameter_name == "topo_discretization")
                    {
                        arret += setTopoDisccretization(line_split[1]);
                    }
                    else if (parameter_name == "transport_diffusion_coefficient")
                    {
                        arret += setTransportDiffusionCoefficient(line_split[1]);
                    }
                    else if (parameter_name == "transport_hydro")
                    {
                        arret += setTransportHydro(line_split[1]);
                    }
                    else if (parameter_name == "v_max")
                    {
                        arret += setVmax(line_split[1]);
                    }
                    //== testing the allowed parameters
                    else
                    {
                        ++arret;
                        logError() << "  - Line" << line_nb
                                   << "Unknown parameter:" << line_split[0];
                        logFatal(__FILE__, __func__, __LINE__);

                    }
                }
            }

            ++line_nb;
        }

        if (arret > 0)
        {
            return EXIT_FAILURE;
        }

        logInfo() << "...done";
    }
    setReferenceList();
    return EXIT_SUCCESS;
}


bool swSolverParam::readParameters()
/**
 * @brief read the parameters files : input_solver_file, extraction_time, probes
 * @return bool of exit_success/failure
 */
{
    logInfo() << "BEGIN IMPORT";
    logInfo() << "Import settings from:" << swFile::input_solver_file;

    if (readInputFile() == EXIT_SUCCESS)
    {
        logInfo() << "Checking parameter consistency";
        if (checkConsistency() == EXIT_SUCCESS)
        {
            logInfo() << "        ...done";
        }
        else
        {
            logError() << "Parameters are inconsistent";
            logFatal(__FILE__, __func__, __LINE__);
            return EXIT_FAILURE;
        }
    }
    else
    {
        logError() << "Error in read()";
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }

    if (readExtractionTime() != EXIT_SUCCESS )
    {
        logError() << "Errors while reading extraction time(s) inout file file";
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }

    logInfo() << "Importing probes info from:" << swFile::input_probe_locations_file;
    if (readProbes() != EXIT_SUCCESS)
    {
        logError() << "Errors while reading probes file";
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    logInfo() << "END IMPORT" << std::endl;
    return EXIT_SUCCESS;
}


bool swSolverParam::readProbes(void)
/**
 * @brief Imports probe values
 * @return bool of exit_success
 */
{
    ///-- set path
    const std::string filename = m_io->absolutePath(swFile::input_probe_locations_file);
    //==  set path

    ///-- Read probe file
    if (std::filesystem::exists(filename))
    {  
        std::ifstream file;
        file.open(filename);
        if (file.is_open() == false)
        {
            logError() << "    Unable to open:" << filename;
            logFatal(__FILE__, __func__, __LINE__);
        }
        else
        {
            size_t line_nb = 0;
            std::string line_txt = swIO::readLine(file);
            ++line_nb;
            line_txt = swIO::readLine(file);
            ++line_nb;
            bool conversion_success = true;
            while (line_txt != "")
            {
                std::vector<std::string> line_split = swString::splitStr(line_txt, swString::regexWhiteSpace());
                if (!(line_split.size() == 3))
                {
                    logError() << " line is misformatted:"  << line_txt;
                    swIO::errorReadingLine(file, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
                    return EXIT_FAILURE;
                }
                else
                {
                    std::string name_probe = line_split[2]; 
                    probe_names.push_back(name_probe);
                    double xprobe = swString::stringToDouble(line_split[0], conversion_success);
                    double yprobe = swString::stringToDouble(line_split[1], conversion_success);
                    if (conversion_success == false)
                    {
                        logError() << "line is misformatted: format conversion";
                        swIO::errorReadingLine(file, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
                    }
                    probe_x.push_back(xprobe);
                    probe_y.push_back(yprobe);
                }
                line_txt = swIO::readLine(file);
                ++line_nb;
            }
            if (file.eof() == false)
            {
                logInfo() << "      Empty line encountered ";
            }
            file.close();
        }
        logInfo() << "      "  << probe_names.size() << " probe(s) read";
    }
    ///==  Read probe file
    return EXIT_SUCCESS;
}


bool swSolverParam::setCFLMax(const std::string& text)
/**
 * @brief Sets the value of the cfl_max (CFL limitation) parameter
 * @param[in] text containing the numerical value for cfl_max
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    cfl_max = swString::stringToDouble(text, conversion_success);//.toDouble(&conversion_success);
    if (!conversion_success ||
        (cfl_max < 0))
    {
        logError() << "  unallowed value for parameter 'cfl_max'" << text;
        logError() << "    'cfl_max' should be positive";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setCoeffRestitution(const std::string& text)
/**
 * @brief Sets the value of the contact_restitution (contact force) parameter
 * @param[in] text containing the numerical value for contact_restitution
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    coeff_restitution = swString::stringToDouble(text, conversion_success);//.toDouble(&conversion_success);
    if (!conversion_success ||
        (coeff_restitution < 0) || (coeff_restitution > 1))
    {
        logError() << "  unallowed value for parameter 'CoeffRestitution'" << text;
        logError() << "    'CoeffRestitution' should be in [0, 1]";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setContactExtMargin(const std::string& text)
/**
 * @brief Sets the value of the contact margin (minimum distance) parameter
 * @param[in] text containing the numerical value for contact_margin
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    contact_ext_margin = swString::stringToDouble(text, conversion_success);//.toDouble(&conversion_success);
    if (!conversion_success ||
        (coeff_restitution < 0))
    {
        logError() << "  unallowed value for parameter 'ContactExtmargin'" << text;
        logError() << "    'CoeffRestitution' should be positif";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setDdpChiPhiMax(const std::string& text)
/**
 * @brief Sets the max value of chi_phi for the DDP model
 * @param[in] text containing the numerical value for ddp_chi_phi_max
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    ddp_chi_phi_max = swString::stringToDouble(text, conversion_success);
    if (ddp_chi_phi_max < 1.0)
    {
        logError() << "  unallowed value for parameter 'ddp_chi_phi_max'" << text;
        logError() << "    allowed values for parameter 'ddp_chi_phi_max' must be larger than unity.";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setDdpChiThetaMax(const std::string& text)
/**
 * @brief Sets the max value of chi_theta for the DDP model
 * @param[in] text containing the numerical value for ddp_chi_theta_max
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    ddp_chi_theta_max = swString::stringToDouble(text, conversion_success);
    if (ddp_chi_theta_max < 1.0)
    {
        logError() << "  unallowed value for parameter 'ddp_chi_theta_max'" << text;
        logError() << "    allowed values for parameter 'ddp_chi_theta_max' must be larger than unity.";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setDdpZtab0Log(const std::string& text)
/**
 * @brief Sets the flag for ztab[0] log file for the DDP model
 * @param[in] text containing the possible options for ddp_ztab0_log
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    int val = swString::stringToInt(text, conversion_success); //.toInt(&conversion_success);
    if (!(conversion_success ||
          (val == 0) ||
          (val == 1)))
    {
        logError() << "  unallowed value for parameter 'ddp_ztab0_log'" << text;
        logError() << "    allowed values for parameter 'ddp_ztab0_log'";
        logError() << "    - 0\t";
        logError() << "    - 1\t";
        return EXIT_FAILURE;
    }
    if (val == 0)
    {
        ddp_ztab0_log = false;
    }
    if (val == 1)
    {
        ddp_ztab0_log = true;
    }

    return EXIT_SUCCESS;
}


bool swSolverParam::setDivergenceCorrection(const std::string& text)
/**
 * @brief Sets the value of the divergence correction
 * @param[in] text containing the method (none, flux limitation, time-step limitation)
 * @return bool of exit_success/failure
 */
{
    divergence_correction = swString::toLower(text);
    if (!((divergence_correction == "none") ||
          (divergence_correction == "limit_flux") ||
          (divergence_correction == "limit_timestep")))
    {
        logError() << "  unallowed value for parameter 'divergence_correction'" << text;
        logError() << "    allowed values for parameter 'divergence_correction'";
        logError() << "    - none\tno divergence correction";
        logError() << "    - limit_flux\tfluxes correction ";
        logError() << "    - limit_timestep\ttimestep correction";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setDtMap(const std::string& text)
/**
 * @brief Sets the value of dt\_map (time period for storing data)
 * @param[in] text containing the numerical value for dt_map
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    dt_map = swString::stringToDouble(text, conversion_success); //.toDouble(&conversion_success);
    if (!conversion_success ||
        (dt_map < -1)) // -1 is authorized to force writing at every time step
    {
        logError() << "  unallowed value for parameter 'dt_map'" << text;
        logError() << "    'dt_map' should be positive or =-1 ";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setDtMax(const std::string& text)
/**
 * @brief Sets the value of the dt_max (maximum value for time-step) prameter
 * @param[in] text containing the numerical value for dt_max
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    dt_max = swString::stringToDouble(text, conversion_success); //.toDouble(&conversion_success);
    if (!conversion_success ||
        (dt_max < 0))
    {
        logError() << "  unallowed value for parameter 'dt_max'" << text;
        logError() << "    'dt_max' should be positive";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setDtTimeSeries(const std::string& text)
/**
 * @brief Sets the value of dt_time_series (time period for storing data for probes)
 * @param[in] text containing the numerical value for dt_time_series
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    dt_time_series = swString::stringToDouble(text, conversion_success);//.toDouble(&conversion_success);
    if (!conversion_success ||
        (dt_time_series < 0))
    {
        logError() << "  unallowed value for parameter 'dt_time_series'" << text;
        logError() << "    'dt_time_series' should be positive";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

bool swSolverParam::setExtractionTimes(std::vector<std::pair<unsigned long, std::string>> &times)
/**
 * @brief Sets the value of extraction_times
 * @param[in] std::vector<std::pair<unsigned long, std::string>> time containing the pairs of time and time string
 * @return bool of exit_success/failure
 */
{
    extraction_times = times;
    return EXIT_SUCCESS;
}


bool swSolverParam::setFrMax(const std::string& text)
/**
 * @brief Sets the value of the Frmax (maximum value for Froude number) prameter
 * @param[in] text containing the numerical value for FrMax
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    Fr_max = swString::stringToDouble(text, conversion_success); //.toDouble(&conversion_success);
    if (!((conversion_success == true) ||
          (Fr_max > 0)))
    {
        logError() << "  unallowed value for parameter 'Fr_max'" << text;
        logError() << "    'Fr_max' should be a positive number";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


bool swSolverParam::setGravity(const std::string& text)
/**
 * @brief Sets the value for the gravity g
 * @param[in] text containing the numerical value for g (must be >0)
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    g = swString::stringToDouble(text, conversion_success);
    if (!conversion_success ||
        (g < 0))
    {
        logError() << "  unallowed value for parameter 'gravity'" << text;
        logError() << "    'gravity' should be positive";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setHmin(const std::string& text)
/**
 * @brief Sets the value of the h_min (minimum non-zero value for water height)
 * @param[in] text containing the numerical value for h_min
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    h_min = swString::stringToDouble(text, conversion_success); // .toDouble(&conversion_success);
    if (!conversion_success ||
        (h_min < 0))
    {
        logError() << "  unallowed value for parameter 'h_min'" << text;
        logError() << "    'h_min' should be positive";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setLogWrite(const std::string& text)
/**
 * @brief Sets the value of the h_min (minimum non-zero value for water height)
 * @param[in] text containing the numerical value for h_min
 * @return bool of exit_success/failure
 */
{
    std::string val = swString::toLower(text);
    if (val == "true")
    {
        log_write = true;
    }
    else if (val == "false")
    {
        log_write = false;
    }
    else
    {
        logError() << "  unallowed value for parameter 'log_write'" << text;
        logError() << "    Allowed values are 'true' or 'false'";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setHotstart(const std::string &text)
/**
 * @brief Sets the value of the hotstart parameter (to restart from a given result file)
 * @param[in] text containing the numerical value for hotstart
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    hotstart_time = swString::stringToDouble(text, conversion_success);
    if (!conversion_success)
    {
        logError() << "  unallowed value for parameter 'hotstart_time'" << text;
        logError() << "    'hotstart_time' should be a number";
        hotstart = false;
        return EXIT_FAILURE;
    }
    hotstart = true;

    return EXIT_SUCCESS;
}


bool swSolverParam::setModelBuilding(const std::string& text)
/**
 * @brief Sets the model for buildings interaction
 * @param[in] text containing the method (none)
 * @return bool of exit_success/failure
 */
{
    /** @todo we don't have any mr-tests on this feature, do we? */
    model_buildings_exchange = swString::toLower(text);
    if (!((model_buildings_exchange == "none") ||
          (model_buildings_exchange == "cell_average")))
    {
        logError() << "  unallowed value for parameter 'model_building'" << text;
        logError() << "    allowed values for parameter 'model_building'";
        logError() << "    'none'\tno building exchange";
        logError() << "    'cell_average'\tcell_average building exchange";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setModelFriction(const std::string& text)
/**
 * @brief Sets the model for friction
 * @param[in] text containing the method (none, Chezy, MAnning, Strickler)
 * @return bool of exit_success/failure
 */
{
    model_friction = swString::toLower(text);
    if (!((model_friction == "none") ||
          (model_friction == "chezy") ||
          (model_friction == "manning") ||
          (model_friction == "strickler")))
    {
        logError() << "  invalid value for parameter 'model_friction'" << text;
        logError() << "  Allowed values are:";
        logError() << "    'none'\t no friction model";
        logError() << "    'chezy'\t Chezy friction model";
        logError() << "    'manning'\t Manning friction model";
        logError() << "    'strickler'\t  Strickler friction model";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


bool swSolverParam::setModelHyperb(const std::string& text)
/**
 * @brief Sets the hyperbolic mode (with or without porosity)
 * @param[in] text containing the model
 * @return bool of exit_success/failure
 */
{
    model_hyp = swString::toLower(text);
    /*
    if ( (model_hyp != "dip") &&
         (model_hyp != "liswes") &&
         (model_hyp != "mp") &&
         (model_hyp != "sp") &&
         (model_hyp != "swes")
        )
    {
        logError() << "  Invalid value for parameter 'model_hyp'" << text;
        logError() << "    allowed values for the 'model_hyp' parameter are";
        logError() << "    - dip\tDual Integral Porosity model";
        logError() << "    - liswes\tLocal Inertial Shallow Water Equations";
        logError() << "    - mp\tMultiple Porosity model";
        logError() << "    - swes\tclassical Shallow Water Equations";
        logError() << "    - sp\tSingle Porosity model";
        return EXIT_FAILURE;
    }
    */
    return EXIT_SUCCESS;
}


bool swSolverParam::setModelInfiltration(const std::string& text)
/**
 * @brief Sets the infiltration model
 * @param[in] text containing the model
 * @return bool of exit_success/failure
 */
{
    model_infiltration = swString::toLower(text);
    if ((model_infiltration != "none") &&
        (model_infiltration != "h_proportional") &&
        (model_infiltration != "horton"))
    {
        /** @todo if sw then ip but porosity = 1 **/
        logError() << "  Invalid value for parameter 'model_infiltration'" << text;
        logError() << "    allowed values for the 'model_infiltration' parameter are";
        logError() << "    - none\t(no infiltration computed)";
        logError() << "    - h_proportional\t(proportional to water depth)";
        logError() << "    - horton\t(Hortonian runoff, Ksat-based)";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setModelLogJam(const std::string& text)
/**
 * @brief Sets the time scheme for the transport part 
 * @param[in] text containing the model
 * @return bool of exit_success/failure
 */
{
    /** @todo we don't have any mr-tests on this feature, do we? */
    model_logjam = swString::toLower(text);
    if (!((model_logjam == "none") ||
          (model_logjam == "rk4") ||
          (model_logjam == "forwardeuler")))
    {
        logError() << "  unallowed value for parameter 'model_logjam'" << text;
        logError() << "    allowed values for parameter 'model_logjam'";
        logError() << "    - none\tno log jam model";
        logError() << "    - ForwardEuler\t log jam model solved with a ForwardEuler scheme";
        logError() << "    - RK4\t log jam model solved with a RK4 scheme";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setModelPorosityRainfall(const std::string& text)
/**
 * @brief Defines how precipitation is to be accounted for in porosity models
 * @param[in] text containing the method (none or active)
 * @return bool of exit_success/failure
 */
{
    model_porosity_rainfall = swString::toLower(text);
    if ( (model_porosity_rainfall != "everywhere") &&
         (model_porosity_rainfall != "water_only") 
        )
    {
        logError() << "  unallowed value for parameter 'model_porosity_rainfall'" << text;
        logError() << "    allowed values are:";
        logError() << "    - everywhere\train falling onto water + land is converted to runoff";
        logError() << "    - water_only\train falling onto land infiltrates instantaneously and is lost";
        logFatal(__FILE__, __func__, __LINE__);
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setModelRunoff(const std::string& text)
/**
 * @brief Sets the model to include precipitation in the simulation
 * @param[in] text containing the method (none or active)
 * @return bool of exit_success/failure
 */
{
    model_runoff = swString::toLower(text);
    if ( (model_runoff != "none") &&
         (model_runoff != "full_rainfall") &&
         (model_runoff != "infiltration_coefficient") &&
         (model_runoff != "runoff_coefficient") &&
         (model_runoff != "scs")
        )
    {
        logError() << "  unallowed value for parameter 'model_runoff'" << text;
        logError() << "    allowed values for parameter 'model_runoff'";
        logError() << "    - none\tno precipitation model";
        logError() << "    - active\tprecipitation activated";
        logFatal(__FILE__, __func__, __LINE__);
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setModelSingularHeadLoss(const std::string& text)
/**
 * @brief Sets the model to include singular head loss in the simulation
 * @param[in] text containing the method (none)
 * @return bool of exit_success/failure
 */
{
    model_singular_head_loss = swString::toLower(text);
    if (!((model_singular_head_loss == "none") ||
          (model_singular_head_loss == "active")))
    {
        logError() << "  unallowed value for parameter 'model_singular_head_loss'" << text;
        logError() << "    allowed values for parameter 'model_singular_head_loss'";
        logError() << "    - none\tno singular head loss model";
        logError() << "    - active\tsingular head loss model";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setModelTransportDiffusion(const std::string& text)
/**
 * @brief Sets the model to include scalar transport in the simulation
 * @param[in] text containing the method (none)
 * @return bool of exit_success/failure
 */
{
    model_transport_diffusion = swString::toLower(text);
    if (!((model_transport_diffusion == "none") ||
        (model_transport_diffusion == "active")))
    {
        logError() << "  unallowed value for parameter 'model_transport_diffusion'" << text;
        logError() << "    allowed values for parameter 'model_transport_diffusion'";
        logError() << "    - none\tno scalar diffusion model";
        logError() << "    - active\tscalar diffusion model";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setModelTransportReaction(const std::string& text)
/**
 * @brief Sets the model to include transport reaction in the simulation
 * @param[in] text containing the method (none or active)
 * @return bool of exit_success/failure
 */
{
    model_transport_reaction = swString::toLower(text);
    if (!(model_transport_reaction == "none" || model_transport_reaction == "active"))
    {
        logError() << "  unallowed value for parameter 'model_transport_reaction'" << text;
        logError() << "    'model_transport_reaction' shoulde be 'none' or 'active'";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setModelWind(const std::string& text)
/**
 * @brief Sets the model to include wind in the simulation
 * @param[in] text containing the method (none, Smith & Banke)
 * @return bool of exit_success/failure
 */
{
    model_wind = swString::toLower(text);
    if (!(model_wind == "none" || model_wind == "smithbanke"))
    {
        logError() << "  unallowed value for parameter 'model_wind'" << text;
        logError() << "    allowed values for parameter 'model_wind'";
        logError() << "    - none\tno wind model";
        logError() << "    - smithbanke\tSmith & Banke wind model";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setMomentumDiffusion(const std::string& text)
/**
 * @brief Sets the model to include momentum diffusion in the simulation
 * @param[in] text containing the method
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    std::string val = swString::toLower(text); 
    if (!(conversion_success ||
        (val == "none") ||
        (val == "velocity_independent") ) )
    {
        logError() << "  unallowed value for parameter 'model_momentum_diffusion'" << text;
        logError() << "    allowed values for parameter 'model_momentum_diffusion'";
        logError() << "    - none\tno momentum diffusion";
        logError() << "    - velocity_independent\tfixed value for all times";
        return EXIT_FAILURE;
    }
    else
    {
        model_momentum_diffusion = val;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setMomentumDiffusionCoefficient(const std::string& text)
/**
 * @brief Sets the model to include momentum diffusion in the simulation
 * @param[in] text containing the method
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    double val = swString::stringToDouble(text, conversion_success);
    if (!conversion_success || val < 0.0)
    {
        logError() << "  unallowed value for parameter 'momentum_diffusion_coefficient'" << text;
        logError() << "    It must be a real, positive number'";
        return EXIT_FAILURE;
    }
    else
    {
        momentum_diffusion_coefficient = val;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setMomentumDissipation(const std::string &text)
/**
 * @brief Sets the model to include momentum dissipation in the simulation
 * @param[in] text containing the method
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    std::string val = swString::toLower(text); 
    if (!(conversion_success &
          ((val == "none") ||
           (val ==  "active"))))
    {
        logError() << "  unallowed value for parameter 'model_momentum_dissipation'" << text;
        logError() << "    allowed values for parameter 'model_momentum_dissipation'";
        logError() << "    - none\tno momentum dissipation";
        logError() << "    - active\tmomentum dissipation computed";
        return EXIT_FAILURE;
    }
    else 
    {
       model_momentum_dissipation = val;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setNegativeDepthDetection(const std::string& text)
/**
 * @brief Sets the variable to decide if a negative depth is searched in the simulation
 * @param[in] text containing the answer yes (1) or no(0)
 * @return bool of exit_success/failure
 */
{
    if (!((text == "none") || (text == "active")))
    {
        logError() << "  unallowed value for parameter 'negative_depth_detection'" << text;
        logError() << "    allowed values for parameter 'negative_depth_detection'";
        logError() << "    - none\tno check";
        logError() << "    - active\tstop the program if the water depth is negative";
        return EXIT_FAILURE;
    }
    negative_depth_detection = text;
    return EXIT_SUCCESS;
}

bool swSolverParam::setNDecimal(const std::string& text)
/**
 * @brief Sets the number of decimals to be displayed in the result maps
 * @param[in] text containing n_decimal
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    n_decimal = swString::stringToInt(text, conversion_success);
    if (!conversion_success ||
        (n_decimal < 0))
    {
        logError() << "  unallowed value for parameter 'n_decimal'" << text;
        logError() << "    'n_decimal' should be a positive integer";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

bool swSolverParam::setNitDivMax(const std::string& text)
/**
 * @brief Sets the maximum iteration number in the divergence correction
 * @param[in] text containing n_it_div_max
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    n_it_div_max = swString::stringToInt(text, conversion_success); // .toInt(&conversion_success);
    if (!conversion_success ||
        (n_it_div_max < 0))
    {
        logError() << "  unallowed value for parameter 'n_it_div_max'" << text;
        logError() << "    'n_it_div_max' should be positive";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setNPorTabMax(const std::string& text)
/**
 * @brief Sets the maximum number of vertical porosity levels in the DDP model
 * @param[in] text containing n_por_tab_max
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    n_por_tab_max = swString::stringToInt(text, conversion_success); // .toInt(&conversion_success);
    if (!conversion_success ||
        (n_por_tab_max < 1))
    {
        logError() << "  unallowed value for parameter 'n_por_tab_max'" << text;
        logError() << "    'n_por_tab_max' should be larger than 1";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setNumScheme(const std::string& text)
/**
 * @brief Sets the numerical scheme used in the simulation
 * @param[in] text containing the method (Godunov, MUSCL-EVR)
 * @return bool of exit_success/failure
 */
{
    num_scheme = swString::toLower(text);
    if ((num_scheme != "godunov") && (num_scheme != "muscl_evr"))
    {
        logError() << "  unallowed value for parameter 'num_scheme'" << text;
        logError() << "    allowed values for parameter 'num_scheme'";
        logError() << "    - none\tno numerical scheme";
        logError() << "    - godunov\tGodunov scheme";
        logError() << "    - muscl_evr\tMUSCL  eigenvector-based linear reconstruction scheme";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setOutputFormat(const std::string &text)
/**
 * @brief Sets the format of output files (binary or ASCII files)
 * @param[in] text containing the value
 * @return bool of exit_success/failure
 */
{
    output_format = swString::toLower(text);
    if (!((output_format == "bin") || (output_format == "txt") || (output_format == "all")))
    {
        logError() << "  unallowed value for parameter 'output' :" << text;
        logError() << "    allowed values for parameter 'output'";
        logError() << "    - bin\tbinary output file generated";
        logError() << "    - txt\ttext output files generated";
        logError() << "    - all\tbinary and text files generated";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

bool swSolverParam::setPhiMin(const std::string &text)
/**
 * @brief Sets the value of the phi_min parameter (minimum phi value for flux/mass computation)
 * @param[in] text containing the numerical value for phi_min
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    phi_min = swString::stringToDouble(text, conversion_success); // .toDouble(&conversion_success);
    if (!conversion_success ||
        (phi_min < 0))
    {
        logError() << "  unallowed value for parameter 'phi_min'" << text;
        logError() << "    'phi_min' should be positive";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

void swSolverParam::setReferenceList(void)
/**
 * @brief Sets up the reference list for citation purposes
 * @param none
 * @return exit success/failure
 */
{
    // -- Leading publication
    biblio_reference_list = {"Caldas et al. (2021) -  https://doi.org/10.1007/978-981-19-1600-7_2 - https://hal.science/hal-03798534v2" };

    //-- Hyperbolic model
    
    if (model_hyp == "swes")
    {
        biblio_reference_list.push_back("SWEs model: A. de Saint-Venant (1871). C.r.hebd.seances Acad.sci., 73 (1871), pp. 148 - 154 (in French) - https://www.biodiversitylibrary.org/item/23729#page/153/mode/1up");
    }
    /*
    if (model_hyp == "ddp")
    {
        biblio_reference_list.push_back("DDP model: Guinot et al. (2018) - https://doi.org.10.1016/j.advwatres.2018.09.014 - https://hal.science/hal-01884110v1");
    }
    */
    else if (model_hyp == "dip")
    {
        biblio_reference_list.push_back("Integral porosity concept: Sanders et al. (2008) - https://doi.org.10.1016/j.jhydrol.2008.08.009");
        biblio_reference_list.push_back("DIP model: Guinot et al. (2017) - https://doi.org.10.1016/j.advwatres.2017.02.009  - https://hal.science/hal-01465071");
    }
    else if (model_hyp == "liswes")
    {
        biblio_reference_list.push_back("LISWEs model: Xia (1994) - https://doi.org/10.1080/00221686.1994.9728358");
    }
    /*
    else if (model_hyp == "mp")
    {
        biblio_reference_list.push_back("MP model: Guinot (2012) - https://doi.org/10.1016/j.advwatres.2011.11.002 - https://hal.science/hal-02269526");
    }
    */
    else if (model_hyp == "sp")
    {
        biblio_reference_list.push_back("SP model: Guinot & Soares-Frazao (2006) - https://doi.org/10.1002/fld.1059");
    }
    //== Hyperbolic model

    //-- Numerical scheme
    if (num_scheme == "godunov")
    {
        biblio_reference_list.push_back("Godunov scheme: Godunov (1959) Matematicheskii Sbornik. Novaya Seriya, 1959, Volume 47(89), Number 3, Pages 271-306 -  https://hal.science/hal-01620642v1");
        biblio_reference_list.push_back("HLL  solver: Harten et al. (1983) - https://doi.org/10.1137/1025002");
        biblio_reference_list.push_back("HLLC solver: Toro et al. (1994) -  https://doi.org/10.1007/BF01414629");
    }
    else if (num_scheme == "muscl_evr")
    {
        biblio_reference_list.push_back("Godunov-type scheme: Godunov (1959) Matematicheskii Sbornik. Novaya Seriya, 1959, Volume 47(89), Number 3, Pages 271-306 -  https://hal.science/hal-01620642v1");
        biblio_reference_list.push_back("MUSCL-EVR scheme: Soares-Frazao & Guinot (2007) - https://doi.org/10.1002/fld.1242 - ");
    }
    //== Numerical scheme

    //-- Runoff model
    if (model_runoff == "scs")
    {
        biblio_reference_list.push_back("SCS model: Mockus (1956) Chapter 10 in U.S.D.A. Soil Conservation Service (SCS) (1956), Hydrology, in National Engineering Handbook, U.S. Dep. of Agric., Washington, D. C");
    }
    //== Runoff model

    //-- Infiltration model
    if (model_infiltration == "horton")
    {
        biblio_reference_list.push_back("Horton (1940) - https://doi.org/10.1029/TR021i002p00541-1");
    }
    //== Infiltration model
}


bool swSolverParam::setReadWallFromFile(const std::string &text)
/**
 * @brief Sets the model to read Wall for transport from a file
 * @param[in] text containing the boolean
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    int val = swString::stringToInt(text, conversion_success); //.toInt(&conversion_success);
    if (!(conversion_success ||
          (val == 0) ||
          (val == 1)))
    {
        logError() << "  unallowed value for parameter 'read_wall_from_file'" << text;
        logError() << "    allowed values for parameter 'read_wall_from_file'";
        logError() << "    - 0\tno file are taking as input";
        logError() << "    - 1\ta file will be used as input";
        return EXIT_FAILURE;
    }
    if (val == 0)
    {
        read_wall_from_file = false;
    }
    if (val == 1)
    {
        read_wall_from_file = true;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setStoreFlag(const std::string& text)
/**
 * @brief Sets the level of storing information in output files
 * @param[in] text containing the level (none, default variables, all variables)
 * @return bool of exit_success/failure
 */
{
    store_flag = swString::toLower(text);
    if (!((store_flag == "none") ||
          (store_flag == "default") ||
          (store_flag == "all")
        ))
    {
        logError() << "  unallowed value for parameter 'StoreFlag'" << text;
        logError() << "    allowed values for parameter 'StoreFlag'";
        logError() << "    - none\tno output txt file";
        logError() << "    - default\tstore only h, q, r";
        logError() << "    - all\tstore all the variables";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setT0(const std::string& text)
/**
 * @brief Sets the inital time for the simulation
 * @param[in] text containing the value for T0
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    t_0 = swString::stringToDouble(text, conversion_success); //.toDouble(&conversion_success);
    if (!conversion_success)
    {
        logError() << "  unallowed value for parameter 't_0'" << text;
        logError() << "    't_0' should be a number";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


bool swSolverParam::setTmax(const std::string& text)
/**
 * @brief Sets the final time for the simulation
 * @param[in] text containing the value for Tmax
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    t_max = swString::stringToDouble(text, conversion_success); // .toDouble(&conversion_success);
    if (!conversion_success)
    {
        logError() << "  unallowed value for parameter 't_max'" << text;
        logError() << "    't_max' should be a number";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


bool swSolverParam::setTopoDisccretization(const std::string& text)
/**
 * @brief sets the type of topography discretization (dicontinuous or blending)
 * @param[in] text containing the type of discretization
 * @return bool of exit_success/failure
 */
{
    topo_discretization = swString::toLower(text);
    if (!((topo_discretization == "discontinuous") ||
          (topo_discretization == "blending")))
    {
        logError() << "  unallowed value for parameter 'topo_discretization'" << text;
        logError() << "    allowed values for parameter 'topo_discretization'";
        logError() << "    - discontinuous\t hllc solver uses only the discontinuous discretized topography";
        logError() << "    - blending\t hllc solver based on mixing between continuous (delta_h) and discontinuous (delta_z) topography";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setTransportDiffusionCoefficient(const std::string& text)
/**
 * @brief Sets the value of the diffusion coefficient por passive scalar transport
 * @param[in] text containing the numerical value
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    transport_diffusion_coefficient = swString::stringToDouble(text, conversion_success); // .toDouble(&conversion_success);
    if (!conversion_success)
    {
        logError() << "  unallowed value for parameter 'transport_diffusion_coefficient'" << text;
        logError() << "    transport_diffusion_coefficient' should be a number";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setTransportHydro(const std::string& text)
/**
 * @brief Sets the value of transport_hydro
 * @param[in] text containing the way to compute (compute, read)
 * @return bool of exit_success/failure
 */
{
    transport_hydro = swString::toLower(text);
    if (!((transport_hydro == "compute") ||
          (transport_hydro == "read")))
    {
        logError() << "  unallowed value for parameter 'transport_hydro'" << text;
        logError() << "    allowed values for parameter 'transport_hydro'";
        logError() << "    - compute\t compute hydro solving the equations";
        logError() << "    - read\t read hydro from a file ";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool swSolverParam::setVmax(const std::string& text)
/**
 * @brief Sets the maximum admissible value for flow velocities
 * @param[in] text containing the value for Vmax
 * @return bool of exit_success/failure
 */
{
    bool conversion_success = true;
    v_max = swString::stringToDouble(text, conversion_success); //.toDouble(&conversion_success);
    if (!((conversion_success == true) ||
          (v_max > 0)))
    {
        logError() << "  unallowed value for parameter 'v_max'" << text;
        logError() << "    'v_max' should be a positive number";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
