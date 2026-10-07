#pragma once

#include <array>
#include <vector>

#include "../../../common/geometry/swMesh.h"
#include "../../common/swSolverParam.h"
#include "../common/swSource.h"
#include "../swParent.h"


class swSwes : public swParent
{
public:
    swSwes(void);
    ~swSwes(void);

    ///-- General methods (public)
    void computeEdgeFactors();
    bool execute() final;
    void importAttributesFromSolverParam(swSolverParam *solver_param) final;
    bool init(swMesh *mesh) final;
    //== General methods


protected:
    ///-- General methods (protected)
    void checkParameterConsistency(void);
    void computeMusclSlopes(const std::vector<double>& cell_x, const std::vector<double>& cell_y, const std::vector<std::vector<size_t>>& cell_neigh_cell_idx,
                               const std::vector<double>& edge_int_midpoint_x, const std::vector<double>& edge_int_midpoint_y,
                               const std::vector<size_t>& edge_int_celL_idx, const std::vector<size_t>& edge_int_celR_idx,
                               const std::vector<double>& cell_averaged_field, std::vector<double>& slope_x, std::vector<double>& slope_y);
    void correctVelocity(void);
    void hypBoundaryFlux(double const& zbR, double const& zR, double const& qR, double const& rR,
                         double const& lambda_plus_R, std::string const& bc_type, double const& bc_val,
                         double& flux1, double& flux2, double& flux3, double& source_right);
    void hypCompute(std::string const& div_corr_option, size_t n_it_div_max, std::string const& num_scheme, double& dt);
    void hypComputeLRStates(std::string const& num_scheme);
    double hypComputeTimeStep(void);
    void hypComputeFluxes(void);
    void hypComputeBalance(double const& dt);
    void hypDivCorr(std::string const& div_corr_option, size_t n_it_div_max, double& dt);
    void hypDivCorrFluxes(size_t n_it_div_max);
    void hypDivCorrTimeStep(void);
    bool hypGodunovComputeLRStates(void);
    double hypHCritical(double const& q);
    void hypHllc(double const& zbL, double const& zbR,
                double const& zL, double const& qL, double const& rL,
                double const& zR, double const& qR, double const& rR,
                double const& lambda_minus_L, double const& lambda_minus_R,
                double const& lambda_plus_L, double const& lambda_plus_R, 
                double& flux1, double& flux2, double& flux3,
                double& source_left, double& source_right);
    void hypHllcBlending(double const& zbL, double const& zbR,
        double const& zL, double const& qL, double const& rL,
        double const& zR, double const& qR, double const& rR,
        double const& lambda_minus_L, double const& lambda_minus_R,
        double const& lambda_plus_L, double const& lambda_plus_R, 
        double& flux1, double& flux2, double& flux3,
        double& source_left, double& source_right);
    void hypHllcImpervious(double const& zbL, double const& zbR,
        double const& zL, double const& qL,
        double const& zR, double const& qR,
        double const& lambda_minus_L, double const& lambda_plus_R, 
        double& flux1, double& flux2, double& flux3,
        double& source_left, double& source_right);
    bool hypMusclEvrComputeLRStates(void);
    double lambdaMinus(double const& u_cell, double const& c);
    double lambdaPlus(double const& u_cell, double const& c);
    void parabComputeExplicit(void);
    void parabComputeScalarExplicit(const std::vector<double>& diffusion_coefficient, const std::vector<double>& water_depth, std::vector<double>& scalar_cell);
    void setFrictionParameters(void);
    void setInfiltrationParameters(void);
    void setInitialConditions(void);
    void setParabParameters(void);
    void setRunoffParameters(void);
    void setSingularHeadLossParameters(void);
    void setTransportDiffusionParameters(void);
    double specificForce(double const& u_cell, double const& h);
    void sourceFriction(std::string const& model_friction, double const& dt);
    void sourceFrictionChezy(double const& dt);
    void sourceFrictionManning(double const& dt);
    void sourceFrictionStrickler(double const& dt);
    void sourceInfiltration(void);
    void sourceRunoff(void);
    void sourceSingularHeadLoss(double const& dt);
    void sourceWind();
    void sourceWindSmithBanke(const std::vector<double>& stationwise_wind_velocity_x_interpolated,
                              const std::vector<double>& stationwise_wind_velocity_y_interpolated);
    void storeResultBoundaryFluxes(double const& t);
    void storeResultExtractionTimes(const std::string& storage_mode, double const& t, const std::string& label);
    void storeResultMaps(const std::string& storage_mode, double const& t);
    void storeResultProbes(double const& t);
    void storeResultQsections(double const& t);
    void transportScalarCompute(void);
    void writeLogEdgeFactors(void);
    void writeLogFriction(void);
    void writeLogInfiltration(void);
    void writeLogMomentumDiffusion(void);
    void writeLogRunoff(void);
    void writeLogSingularHeadLoss(void);
    void writeLogSolver(void);
    void writeLogWind(void);
    void writeLogs(void);
    //== General methods


    ///-- Protected attributes imported from solver_param
    double cfl_max;
    std::string div_corr_option;
    double dt_max;
    std::vector<std::pair<unsigned long, std::string>>  extraction_times;
    double fr_max;
    bool has_qsection;
    double h_min;
    bool log_write;
    std::string map_storage_mode;
    std::string model_friction;
    std::string model_hyp;
    std::string model_infiltration;
    std::string model_momentum_diffusion;
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
    bool do_extraction = false;
    size_t extraction_time_idx = 0;
    std::vector<size_t> probe_cell_idx = {};
    size_t probes_nb = 0;
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


    //-- Conserved variables
    std::vector<double> scs_h_a = {};                  /**< [i]: depth in the initial accumulation compartment (SCS model) for cell i */
    std::vector<double> scs_h_s = {};                  /**< [i]: depth in the storage compartment (SCS model) for cell i */
    std::vector<double> unit_discharge_x = {};         /**< global x unit-discharge vector (to be used cellwise)*/
    std::vector<double> unit_discharge_x_old = {};     /**< copy of unit_discharge_x before balance (for divergence correction) */
    std::vector<double> unit_discharge_y = {};         /**< global y unit-discharge vector (to be used cellwise)*/
    std::vector<double> unit_discharge_y_old = {};     /**< copy of unit_discharge_y before balance (for divergence correction) */
    std::vector<double> water_depth = {};              /**< global water depth vector (to be used cellwise)*/
    std::vector<double> water_depth_old = {};          /**< copy of water depth before balance (for divergence correction) */
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
    std::vector<double> h_riemann_bnd = {};                     /**< [i]: water depth (Riemann state) on RHS of i-th boundary edge */
    std::vector<std::array<double, 2>> h_riemann_int = {};      /**< [i][0/1]: water depth (Riemann state) on LHS/RHS of i-th internal edge */
    std::vector<double> lambda_max_bnd = {};                    /**< [i]: max. lambda value on RHS of i-th boundary edge */
    std::vector<double> lambda_max_int = {};                    /**< [i]: max. lambda value on both sides of i-th internal edge */
    std::vector<double> lambda_minus_bnd = {};                  /**< [i]: lambda_minus value on RHS of i-th boundary edge */
    std::vector<std::array<double, 2>> lambda_minus_int = {};   /**< [i][0/1] lambda_minus value on LHS/RHS of i-th internal edge */
    std::vector<double> lambda_plus_bnd = {};                   /**< [i]: lambda_plus value on RHS of i-th boundary edge */
    std::vector<std::array<double, 2>> lambda_plus_int = {};    /**< [i][0/1] lambda_plus value on LHS/RHS of i-th internal edge */
    std::vector<double> q_riemann_bnd = {};                     /**< [i]: normal unit discharge (Riemann state) on RHS of i-th boundary edge */
    std::vector<std::array<double, 2>> q_riemann_int = {};      /**< [i][0/1]: normal unit discharge (Riemann state) on LHS/RHS of i-th internal edge */
    std::vector<double> r_riemann_bnd = {};                     /**< [i]: tangent unit discharge (Riemann state) on RHS of i-th boundary edge */
    std::vector<std::array<double, 2>> r_riemann_int = {};      /**< [i][0/1]: tangent unit discharge (Riemann state) on LHS/RHS of i-th internal edge */
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


    ///-- Optional model parameters
    double momentum_diffusion_coefficient = 1e-6;                     /**< Momentum diffusion coefficient */
    std::vector<double> edge_int_diff_coeff_normal = {};              /**< [i] Diffusion coefficient in the direction normal to the i-th interface */
    std::vector<std::vector<size_t>> diffusion_matrix_columns = {};   /**< [i][j] Index of the j-th non-zero column in the i-th row of the diffusion matrix */
    std::vector<std::vector<double>> diffusion_matrix_values = {};    /**< [i][j] Numerical value of the j-th non-zero element in the i-th row of the diffusion matrix */

    std::vector<double> friction_alpha = {};                          /**< [i] friction angle in cell i */
    std::vector<double> friction_coeff_direction_1 = {};              /**< [i] friction coeff in x-direction for cell i */
    std::vector<double> friction_coeff_direction_2 = {};              /**< [i] friction coeff in y-direction for cell i */
    std::vector<double> friction_cos_alpha = {};                      /**< [i] friction cos(angle) in cell i */
    std::vector<double> friction_sin_alpha = {};                      /**< [i] friction sin(angle) in cell i */
    
    std::vector<double> infiltration_param = {};                      /**< [i] value for the infiltration parameter in cell i */
    
    std::vector<double> runoff_param = {};                            /**< [i] value for the runoff parameter in cell i (if model != SCS) */
    std::vector<double> scs_i_a = {};                                 /**< [i] value for the SCS' Initial accumulation Ia in cell i (if model == SCS) */
    std::vector<double> scs_s = {};                                   /**< [i] value for the SCS' Storage capacity S in cell i (if model == SCS) */


    std::vector<double> singular_head_loss_alpha = {};                /**< [i]: angle of the 1st principal direction of singular loss tensor with x-direction in cell i */
    std::vector<double> singular_head_loss_h_sing = {};               /**< [i]: threshold value for sinh. h. loss in cell i */
    std::vector<double> singular_head_loss_coeff_direction_1 = {};    /**< [i]: sing head loss coeff for 1st principal direction in cell i */
    std::vector<double> singular_head_loss_coeff_direction_2 = {};    /**< [i]: sing head loss coeff for 2nd principal direction in cell i */

    std::vector<double> transport_scalar_degradation_parameter = {};    /**< [i]: degradation constant in cell i */
    std::vector<double> transport_scalar_diffusion_parameter = {};    /**< [i]: diffusion or dispersion coefficient in cell i */
    //== Optional model parameters
};
