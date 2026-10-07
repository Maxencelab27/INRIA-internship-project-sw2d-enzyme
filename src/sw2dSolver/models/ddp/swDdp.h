#pragma once

#include <array>
#include <vector>

#include "../../../common/geometry/swMesh.h"
#include "../../common/swSolverParam.h"
#include "../common/swSource.h"
#include "../swParent.h"


class swDdp : public swParent
{
public:
    swDdp(void);
    ~swDdp(void);

    ///-- General methods (public)
    void computeEdgeFactors();
    bool execute() final;
    void importAttributesFromSolverParam(swSolverParam *solver_param) final;
    bool init(swMesh *mesh) final;
    //== General methods


protected:
    ///-- General methods (protected)
    void checkParameterConsistency(void);
    bool computeMusclSlopes(const std::vector<double>& cell_x, const std::vector<double>& cell_y, const std::vector<std::vector<size_t>>& cell_neigh_cell_idx,
                               const std::vector<double>& edge_int_midpoint_x, const std::vector<double>& edge_int_midpoint_y,
                               const std::vector<size_t>& edge_int_celL_idx, const std::vector<size_t>& edge_int_celR_idx,
                               const std::vector<double>& cell_averaged_field, std::vector<double>& slope_x, std::vector<double>& slope_y);   
    bool correctVelocity(void);
    bool hypBoundaryFlux(const size_t& idx_edge_bnd, const double& bR, const double& zbR, const double& zR, const double& qR, const double& rR,
                         const double& lambda_plus_R, std::string const& bc_type, const double& bc_val,
                         double& flux1, double& flux2, double& flux3, double& source_right);
    bool hypCompute(std::string const& div_corr_option, size_t n_it_div_max, std::string const& num_scheme, double& dt);
    bool hypComputeLRStates(std::string const& num_scheme);
    bool hypComputeFluxes(void);
    bool hypComputeBalance(const double& dt);
    bool hypDivCorr(std::string const& div_corr_option, size_t n_it_div_max, double& dt);
    bool hypDivCorrFluxes(size_t n_it_div_max);
    bool hypDivCorrTimeStep(void);
    bool hypGodunovComputeLRStates(void);
    double hypHCriticalPhiq(double const& q, double const& PhiG, double const& PhiW);
    double hypHCriticalPhiWq(double const& phiwq, double const& phiG);
    bool hypHllc(const size_t& idx_edge, double const& bL, double const& bR,
                         double const& zL, double const& uL, double const& vL,
                         double const& zR, double const& uR, double const& vR, 
                         double& flux1, double& flux2, double& flux3,
                         double& source_left, double& source_right);
    bool hypMusclEvrComputeLRStates(void);
    double hypTimeStep(void);
    double lambdaMinus(double const& u_cell, double const& c, double const& x_phi, double const& x_theta, double const& b);
    double lambdaPlus(double const& u_cell, double const& c, double const& x_phi, double const& x_theta, double const& b);
    void parabComputeExplicit(void);
    void parabComputeScalarExplicit(const std::vector<double>& diffusion_coefficient, const std::vector<double>& water_depth, std::vector<double>& scalar_cell);
    double phiGBndOfZ(const size_t& idx_edge_bnd, const double& z);
    double phiGIntOfZ(const size_t& idx_edge_int, const double& z);
    double phiWOfZ(const size_t& idx_cell, const double& z);
    size_t porosityLineSize(const size_t& line_nb, const double& porosity_law_type, const std::vector<std::string> line_txt);
    double pressureForceGammaOfZ(const size_t& idx_edge, const double& z);
    double pressureForceOmegaOfZ(const size_t& idx_cell, const double& z);
    bool computePressureForceOmegaOntoGammaOfZ(const size_t& idx_edge, const double& zL, const double& zR, double& Fp_L, double& Fp_R, double& Fb_L, double& Fb_R);
    double pressureWaveCelerity(double const &theta, double const &phi);
    bool readPorosityFile(const std::string& file_name, const size_t& cell_nb, std::vector<std::string>& porosity_headers, std::vector<std::vector<double>>& porosity_values);
    void setBoussinesqParameter(void);
    void setFrictionParameters(void);
    void setInfiltrationParameters(void);
    void setInitialConditions(void);
    void setMomentumDissipationParameters(void);
    void setParabParameters(void);
    bool setPorosityCellMap(const std::vector<std::vector<double>>& records);
    bool setPorosityEdgeBndMap(const std::vector<std::vector<double>>& records);
    bool setPorosityEdgeIntMap(const std::vector<std::vector<double>>& records);
    bool setPorosityParameters(void);
    void setRunoffParameters(void);
    void setSingularHeadLossParameters(void);
    void setTransportDiffusionParameters(void);
    double specificForce(const double &theta_velocity, const double &velocity, const double &pressure_force);
    void sourceFriction(std::string const& model_friction, const double& dt);
    void sourceFrictionChezy(const double& dt);
    void sourceFrictionManning(const double& dt);
    void sourceFrictionStrickler(const double& dt);
    void sourceInfiltration(void);
    void sourceRunoffEverywhere(void);
    void sourceRunoffWateronly(void);
    void sourceSingularHeadLoss(const double& dt);
    void sourceWind();
    void sourceWindSmithBanke(const std::vector<double>& stationwise_wind_velocity_x_interpolated,
                              const std::vector<double>& stationwise_wind_velocity_y_interpolated);
    void storeResultMaps(const std::string& storage_mode, const double& t);
    void storeResultProbes(const double& t);
    bool tabulatePorosity(size_t ntab, const std::vector<double> &z, const std::vector<double> &phi, 
                                    std::vector<double>& z_tab, std::vector<double>& phi_tab);
    double thetaGBndOfZ(const size_t& idx_edge_bnd, const double& z);
    double thetaGIntOfZ(const size_t& idx_edge_int, const double& z);
    double thetaWOfZ(const size_t& idx_cell, const double& z);
    void writeLogBoussinesq(void);
    void writeLogEdgeFactors(void);
    void writeLogFriction(void);
    void writeLogInfiltration(void);
    void writeLogMomentumDiffusion(void);
    bool writeLogPorosity(void);
    void writeLogRunoff(void);
    void writeLogSingularHeadLoss(void);
    void writeLogSolver(void);
    void writeLogWind(void);
    void writeLogs(void);
    double zOfThetaW(const size_t& idx_cell, const double& theta);
    //== General methods


    ///-- Protected attributes imported from solver_param
    double cfl_max;
    double chi_phi_max;
    double chi_theta_max;
    std::string div_corr_option;
    double dt_map;
    double dt_max;
    double dt_time_series;
    double fr_max;
    double h_min;
    double h_min_ddp;                                         
    bool log_write;
    std::string map_storage_mode;
    std::string model_friction;
    std::string model_hyp;
    std::string model_infiltration;
    std::string model_momentum_diffusion;
    std::string model_momentum_dissipation;
    std::string model_porosity_rainfall;
    std::string model_runoff;
    std::string model_singular_head_loss;
    std::string model_transport_diffusion;
    std::string model_wind;
    size_t momentum_diff_max_iter;
    double momentum_diff_convergence_eps;
    std::string negative_depth_detection;
    std::string num_scheme;
    size_t n_decimal;
    size_t n_it_div_max;
    size_t n_por_tab_max;
    std::string output_format;
    std::vector<std::string> probe_names;
    std::vector<double> probe_x;
    std::vector<double> probe_y;
    double t_0;
    double t_max;
    std::string topo_discretization;
    double transport_diffusion_coefficient;
    double v_max;
    //== Protected attributes imported from solver_param

    ///-- Protected local Attributes
    double dt = 0.;
    std::vector<size_t> probe_cell_idx = {};
    size_t probes_nb = 0;
    bool store_map = false;
    bool store_time_series = false;
    double t = 0.0;
    double t_cpu_beg = 0.;
    double total_t_clock_hyperb = 0.;
    double total_t_clock_hyperb_lrstates = 0.;
    double total_t_clock_hyperb_fluxes = 0.;
    double total_t_clock_hyperb_balance = 0.;
    double total_t_clock_hyperb_divcorr = 0.;
    double total_t_clock_hyperb_timestep = 0.;
    double total_t_clock_parab = 0.;
    double total_t_clock_source = 0.;
    double total_t_clock_storage = 0.;
    size_t warning_thetaG_gt_thetaW = 0;
    //== Protected local Attributes


    ///-- Conserved variables
    std::vector<double> scs_h_a = {};                  /**< [i]: depth in the initial accumulation compartment (SCS model) for cell i */
    std::vector<double> scs_h_s = {};                  /**< [i]: depth in the storage compartment (SCS model) for cell i */
    std::vector<double> thetaG_bnd = {};               /**< [i] state variable thetaG at i-th boundary edge*/
    std::vector<double> thetaG_int = {};               /**< [i] state variable thetaG at i-th internal edge*/
    std::vector<double> thetaW = {};                   /**< [i] state variable thetaW at i-th cell*/
    std::vector<double> thetaW_velocity_x = {};        /**< [i] state variable thetaW*u at i-th cell*/
    std::vector<double> thetaW_velocity_y = {};        /**< [i] state variable thetaW*v at i-th cell*/
    //== Conserved variables


    ///-- Fluxes & interface-based source terms
    std::vector<double> F1_bnd = {};                   /**< [i]: normal mass flux (hu) across i-th boundary edge */
    std::vector<double> F2_bnd = {};                   /**< [i]: normal momentum flux across i-th boundary edge */
    std::vector<double> F3_bnd = {};                   /**< [i]: tangent momentum flux across i-th boundary edge */
    std::vector<double> F1_int = {};                   /**< [i]: normal mass flux (hu) across i-th internal edge */
    std::vector<double> F2_int = {};                   /**< [i]: normal momentum flux across i-th internal edge */
    std::vector<double> F3_int = {};                   /**< [i]: tangent momentum flux across i-th internal edge */
    std::vector<double> s2L_int = {};                  /**< [i]: momentum source term on LHS of i-th internal interface */
    std::vector<double> s2R_bnd = {};                  /**< [i]: momentum source term on RHS of i-th boundary interface */
    std::vector<double> s2R_int = {};                  /**< [i]: momentum source term on RHS of i-th internal interface */
    //== Fluxes & interface-based source terms


    ///-- Interface-based Riemann states
    std::vector<double> lambda_max_bnd = {};                    /**< [i]: max. lambda value on RHS of i-th boundary edge */
    std::vector<double> lambda_max_int = {};                    /**< [i]: max. lambda value on both sides of i-th internal edge */
    std::vector<double> lambda_minus_bnd = {};                  /**< [i]: lambda_minus value on RHS of i-th boundary edge */
    std::vector<std::array<double, 2>> lambda_minus_int = {};   /**< [i][0/1] lambda_minus value on LHS/RHS of i-th internal edge */
    std::vector<double> lambda_plus_bnd = {};                   /**< [i]: lambda_plus value on RHS of i-th boundary edge */
    std::vector<std::array<double, 2>> lambda_plus_int = {};    /**< [i][0/1] lambda_plus value on LHS/RHS of i-th internal edge */
    std::vector<double> theta_riemann_bnd = {};                 /**< [i]: water depth (Riemann state) on RHS of i-th boundary edge */
    std::vector<std::array<double, 2>> theta_riemann_int = {};  /**< [i][0/1]: water depth (Riemann state) on LHS/RHS of i-th internal edge */
    std::vector<double> theta_u_riemann_bnd = {};               /**< [i]: normal unit discharge (Riemann state) on RHS of i-th boundary edge */
    std::vector<std::array<double, 2>> theta_u_riemann_int = {};/**< [i][0/1]: normal unit discharge (Riemann state) on LHS/RHS of i-th internal edge */
    std::vector<double> theta_v_riemann_bnd = {};               /**< [i]: tangent unit discharge (Riemann state) on RHS of i-th boundary edge */
    std::vector<std::array<double, 2>> theta_v_riemann_int = {};/**< [i][0/1]: tangent unit discharge (Riemann state) on LHS/RHS of i-th internal edge */
    std::vector<double> z_riemann_bnd = {};                     /**< [i]: surface elevation (Riemann state) on RHS of i-th boundary edge */
    std::vector<std::array<double, 2>> z_riemann_int = {};      /**< [i][0/1]: surface elevation (Riemann state) on LHS/RHS of i-th internal edge */
    //== Interface-based Riemann states


    ///-- Edge factors for faster balance
    std::vector<double> edge_bnd_factor_1_F1 = {};     /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_bnd_factor_2_F2 = {};     /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_bnd_factor_2_F3 = {};     /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_bnd_factor_2_S2 = {};     /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_bnd_factor_3_F2 = {};     /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_bnd_factor_3_F3 = {};     /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_bnd_factor_3_S2 = {};     /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_1_F1_L = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_1_F1_R = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_2_F2_L = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_2_F2_R = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_2_F3_L = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_2_F3_R = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_2_S2_L = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_2_S2_R = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_3_F2_L = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_3_F2_R = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_3_F3_L = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_3_F3_R = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_3_S2_L = {};   /**< pre-compute factors for a faster computation of fluxes*/
    std::vector<double> edge_int_factor_3_S2_R = {};   /**< pre-compute factors for a faster computation of fluxes*/
    //== Edge factors for faster balance


    ///-- Model parameters
    std::vector<double> boussinesq = {};                             /**< [i] Boussinesq coefficient in the i-th cell */
    std::vector<double> cell_mu_1 = {};                               /**< [i] momentum dissipation coefficient for 1st principal direction in cell i */
    std::vector<double> cell_mu_2 = {};                               /**< [i] momentum dissipation coefficient for 2nd principal direction in cell i */
    std::vector<double> cell_mu_alpha = {};                           /**< [i] angle of 1st principal direction with x-axis in cell i */
    std::vector<double> cell_mu_cos_alpha = {};                       /**< [i] cos(angle) of 1st principal direction with x-axis in cell i */
    std::vector<double> cell_mu_sin_alpha = {};                       /**< [i] sin(angle) of 1st principal direction with x-axis in cell i */
    std::vector<std::vector<double>> cell_z_tab = {};                 /**< [i][k] vertical level k for in the i-th cell porosity */
 
    double momentum_diffusion_coefficient = 1e-6;                     /**< Momentum diffusion coefficient */

    std::vector<std::vector<double>> edge_bnd_z_tab = {};             /**< [i][k] vertical level k for in the i-th boundary edge porosity */
    std::vector<double> edge_int_diff_coeff_normal = {};              /**< [i] Diffusion coefficient in the direction normal to the i-th interface */
    std::vector<std::vector<double>> edge_int_z_tab = {};             /**< [i][k] vertical level k for in the i-th internal edge porosity */
    std::vector<std::vector<size_t>> diffusion_matrix_columns = {};   /**< [i][j] Index of the j-th non-zero column in the i-th row of the diffusion matrix */
    std::vector<std::vector<double>> diffusion_matrix_values = {};    /**< [i][j] Numerical value of the j-th non-zero element in the i-th row of the diffusion matrix */

    std::vector<double> friction_alpha = {};                          /**< [i] friction angle in cell i */
    std::vector<double> friction_coeff_direction_1 = {};              /**< [i] friction coeff in x-direction for cell i */
    std::vector<double> friction_coeff_direction_2 = {};              /**< [i] friction coeff in y-direction for cell i */
    std::vector<double> friction_cos_alpha = {};                      /**< [i] friction cos(angle) in cell i */
    std::vector<double> friction_sin_alpha = {};                      /**< [i] friction sin(angle) in cell i */
    
    std::vector<double> infiltration_param = {};                      /**< [i] value for the infiltration parameter in cell i */
    std::vector<std::vector<double>> phiG_bnd_tab = {};               /**< [i][k] connectivity porosity over the i-th boundary edge, at vertical level k */
    std::vector<std::vector<double>> phiG_int_tab = {};               /**< [i][k] connectivity porosity over the i-th internal edge, at vertical level k */
    std::vector<std::vector<double>> phiW_tab = {};                   /**< [i][k] storage porosity in the i-th cell, at vertical level k */
    
    std::vector<double> runoff_param = {};                            /**< [i] value for the runoff parameter in cell i (if model != SCS) */
    std::vector<double> scs_i_a = {};                                 /**< [i] value for the SCS' Initial accumulation Ia in cell i (if model == SCS) */
    std::vector<double> scs_s = {};                                   /**< [i] value for the SCS' Storage capacity S in cell i (if model == SCS) */

    std::vector<double> singular_head_loss_alpha = {};                /**< [i]: angle of the 1st principal direction of singular loss tensor with x-direction in cell i */
    std::vector<double> singular_head_loss_h_sing = {};               /**< [i]: threshold value for sinh. h. loss in cell i */
    std::vector<double> singular_head_loss_coeff_direction_1 = {};    /**< [i]: sing head loss coeff for 1st principal direction in cell i */
    std::vector<double> singular_head_loss_coeff_direction_2 = {};    /**< [i]: sing head loss coeff for 2nd principal direction in cell i */

    std::vector<std::vector<double>> thetaG_bnd_tab = {};             /**< [i][k] connectivity integral porosity over the i-th boundary edge, at vertical level k */
    std::vector<std::vector<double>> thetaG_int_tab = {};             /**< [i][k] connectivity integral porosity over the i-th internal edge, at vertical level k */
    std::vector<double> thetaG_int_hmin = {};                         /**< [i] minimal value of theta_gamma in the HLLC computation */
    std::vector<double> thetaW_old = {};                              /**< copy of state variable thetaW before balance (for divergence correction & momentum dissipation) */
    std::vector<std::vector<double>> thetaW_tab = {};                 /**< [i][k] storage integral porosity in the i-th cell, at vertical level k */
    std::vector<double> thetaW_velocity_x_old = {};                   /**< copy of state variable thetaW*u before balance (for divergence correction) */
    std::vector<double> thetaW_velocity_y_old = {};                   /**< copy of state variable thetaW*v before balance (for divergence correction) */
    //== Model parameters
};
