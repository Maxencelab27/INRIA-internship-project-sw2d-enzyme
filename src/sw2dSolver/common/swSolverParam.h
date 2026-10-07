#pragma once

#include <string>
#include <vector>

#include "../../common/param.h"

class swIO;

class swSolverParam
{
public:
    swSolverParam(swIO *io);
    ~swSolverParam() = default;

    bool checkConsistency();
    bool readExtractionTime();
    bool readInputFile();
    bool readParameters();
    bool readProbes(void);

    ///-- Getters
    double getCFLMax();
    std::vector<std::string> getBiblioReferenceList(); 
    std::string getDivergenceCorrection();
    double getDtMap();
    double getDtMax();     
    double getDtTimeSeries();
    std::vector<std::pair<unsigned long, std::string>>& getExtractionTimes();
    double getFrMax();
    double getGravity();
    double getHMin();
    double getHMinDdp();
    bool getHasQsection();
    bool getHotstart();
    bool getLogWrite();
    std::string getModelBuildingsExchange();
    std::string getModelFriction();
    std::string getModelHyp();
    std::string getModelInfiltration();
    std::string getModelLogJam();
    std::string getModelMomentumDiffusion();
    std::string getModelMomentumDissipation();
    std::string getModelParab();
    std::string getModelPorosityRainfall();
    std::string getModelRunoff();
    std::string getModelSingularHeadLoss();
    std::string getModelTransportReaction();
    std::string getModelWind();
    double getMomentumDiffusionCoefficient();
    size_t getMomentumDiffMaxIter();
    double getMomentumDiffConvergenceEps();
    bool getMomentumDissipation();
    size_t getNDecimal();
    size_t getNItDivMax();
    size_t getNPorTabMax();
    std::string getNegativeDepthDetection();
    std::string getNumScheme();
    std::string getOutputFormat();
    double getPhiMin();
    std::vector<std::string> getProbeNames();
    std::vector<double> getProbeX();
    std::vector<double> getProbeY();
    std::string getStoreFlag();
    double getT0();
    double getTMax();
    std::string getTopoDiscretization();
    double getTransportDiffusionCoefficient();
    double getVMax();
    //=== Getters

    ///-- Setters
    bool setCFLMax(const std::string& text);
    double getDdpChiPhiMax();
    double getDdpChiThetaMax();
    std::string getDdpLawType();
    bool setCoeffRestitution(const std::string& text);
    bool setContactExtMargin(const std::string& text);
    bool setDdpChiPhiMax(const std::string& text);
    bool setDdpChiThetaMax(const std::string& text);
    bool setDdpZtab0Log(const std::string& text);
    bool setDivergenceCorrection(const std::string& text);
    bool setDtMap(const std::string& text);
    bool setDtMax(const std::string& text);
    bool setDtTimeSeries(const std::string& text);
    bool setExtractionTimes(std::vector<std::pair<unsigned long, std::string>> &times);
    bool setFrMax(const std::string& text);
    bool setGravity(const std::string& text);
    bool setHasQsection(const std::string &text);
    bool setHmin(const std::string& text);
    bool setHotstart(const std::string& text);
    bool setLogWrite(const std::string& text);
    bool setModelBuilding(const std::string& text);
    bool setModelFriction(const std::string& text);
    bool setModelInfiltration(const std::string& text);
    bool setModelGates(const std::string& text);
    bool setModelHyperb(const std::string& text);
    bool setModelLogJam(const std::string& text);
    bool setModelOrifice(const std::string& text);
    bool setModelPorosityRainfall(const std::string& text);
    bool setModelRunoff(const std::string& text);
    bool setModelPump(const std::string& text);
    bool setModelSingularHeadLoss(const std::string& text);
    bool setModelTransportDiffusion(const std::string& text);
    bool setModelTransportReaction(const std::string& text);
    bool setModelWind(const std::string& text);
    bool setMomentumDiffusion(const std::string& text);
    bool setMomentumDiffusionCoefficient(const std::string& text);
    bool setMomentumDissipation(const std::string& text);
    bool setNegativeDepthDetection(const std::string& text);
    bool setNDecimal(const std::string& text);
    bool setNitDivMax(const std::string& text);
    bool setNPorTabMax(const std::string& text);
    bool setNumScheme(const std::string& text);
    bool setOutputFormat(const std::string& text);
    bool setPhiMin(const std::string& text);
    void setReferenceList(void);
    bool setReadWallFromFile(const std::string& text);
    bool setStoreFlag(const std::string& text);
    bool setT0(const std::string& text);
    bool setTmax(const std::string& text);
    bool setTopoDisccretization(const std::string& text);
    bool setTransportDiffusionCoefficient(const std::string& text);
    bool setTransportHydro(const std::string& text);
    bool setVmax(const std::string& text);
    //--- Setters

public:
    ///-- Public Attributes
    swIO *m_io = nullptr;                   /**< pointer to the swIO object to handle files IO */

private:
    ///-- Private Attributes
    // Extraction times
    std::vector< std::pair<unsigned long,std::string>> extraction_times;

    //// model names for the various processes
    std::string model_buildings_exchange = "none";
    std::string model_friction = "none";
    std::string model_hyp = "swes";
    std::string model_infiltration = "none";
    std::string model_logjam = "none";
    std::string model_momentum_diffusion = "none";
    std::string model_momentum_dissipation = "none";
    std::string model_parab = "none";
    std::string model_porosity_rainfall = "water_only";
    std::string model_runoff = "none";
    std::string model_singular_head_loss = "none";
    std::string model_transport_diffusion = "none";
    std::string model_transport_reaction = "none";
    std::string model_wind = "none";
    std::string num_scheme = "godunov";             /**< Numerical scheme used for the hyperbolic part **/
    std::string topo_discretization = "discontinuous";


    //// general attributes
    double cfl_max = 1;                             /**< max. permissible CFL**/
    std::vector<std::string> biblio_reference_list; /**< [i] The i-th bibliographical reference in the list. */
    double contact_ext_margin = 0.01;               /**< exterior contact for margin**/
    double coeff_restitution = 1;                   /**< coefficient of restitution (contact force) (previously contact_restitution) **/
    double ddp_chi_phi_max = 20.0;                  /**< maximum value for \xi_phi in the DDP model**/
    double ddp_chi_theta_max = 20.0;                /**< maximum value for \xi_theta in the DDP model**/
    bool ddp_ztab0_log = false;                     /**< @todo should the ztab log file be written?**/
    std::string divergence_correction = "none";     /**< divergence correction flag**/
    double dt_map = 100;                            /**< map result storage period**/
    double dt_max = 10;                             /**< max. permissible time step**/
    double dt_time_series = 100;                    /**< time series storage period**/
    double Fr_max = 100;                            /**< maximal allowed Froude**/
    double g = 9.81;                                /**< gravitational acceleration**/
    double h_min = 1e-5;                            /**< water depth for wet/dry criteria**/
    double h_min_ddp = 1e-3;                        /**< h minimal value to compute thetaG_int_hmin (used in HLLC for DDP) */
    bool has_qsection = false;                      /**< flag to compute q sections  */
    bool hotstart = false;                          /**< flag for hotstart**/
    double hotstart_time;                           /**< time for hotstart**/
    bool log_write = false;                         /**< flag for log file writing (= false saves time and disk space in case of big models) */
    double momentum_diffusion_coefficient = 1e-6;   /**< Momentum diffusuon coefficient (by default: molecular diffusion in water */
    double momentum_diff_convergence_eps = 1e-6;    /**< Precision (residual) for conjugate gradient method **/
    size_t momentum_diff_max_iter = 10;             /**< Maximum iteration number for momentum diffusion implicit scheme **/
    size_t n_decimal = 10;                          /**< number of decimals to display in the result maps */
    size_t n_it_div_max = 10;                       /**< max. no. iterations in divergence correction**/
    size_t n_por_tab_max = 5;                       /**< max. no. of vertical porosity levels in DDP model**/
    std::string negative_depth_detection = "none";  /**< flag to stop simulation if negative water depth is detected*/
    std::string output_format = "all";              /**< specifies the output format **/
    double phi_min = small_real;                    /**< phi threshold for flux/balance computation **/
    std::vector<std::string> probe_names = {};      /**< vector containing all probe names **/
    std::vector<double> probe_x = {};               /**< vector containing x-position of every probe**/
    std::vector<double> probe_y = {};               /**< vector containing y-position of every probe**/
    bool read_wall_from_file = false;
    std::string store_flag = "default";             /**< specifies the level of detail of result storage**/
    double t_0 = 0;                                 /**< initial simu. time**/
    double t_max = 100;                             /**< final simu. time**/
    std::string transport_hydro = "compute";        /**< hydrodynamics is read or computed for passive transport **/
    double transport_diffusion_coefficient = 1e-6;  /**< Transport diffusion coefficient for passive scalar transport */
    double v_max = 100;                             /**< maximal allowed velocity norm**/
    //== Attributes
};
