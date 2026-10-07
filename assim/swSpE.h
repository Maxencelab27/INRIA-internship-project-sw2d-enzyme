#pragma once
// =====================================================================
//  swSpE — couche d'assimilation / différentiation automatique (Enzyme)
//  DÉCOUPLÉE du SW2D natif.
//
//  swSpE dérive de la classe native swSp SANS modifier aucun fichier natif.
//  Possible car swSp n'a aucune section private: (public:/protected:) et
//  swParent est tout public: -> la classe dérivée voit tout l'état natif.
//
//  RÈGLE C++ : une méthode de swSpE (même static) n'accède aux membres
//  protected du natif qu'à travers un swSpE* (pas un swSp*). Les fonctions
//  statiques différenciées par Enzyme prennent donc swSpE*.
//
//  Le natif reste STRICTEMENT pristine ; toutes les régularisations
//  Enzyme-safe vivent ici.
// =====================================================================

#include <array>
#include <vector>
#include <map>
#include <string>
#include <cmath>

#include "../src/sw2dSolver/models/sp/swSp.h"
#include "../src/common/param.h"
#include "../src/common/swBoundaryCondition.h"                       // type de base (getBCType)
#include "../src/sw2dSolver/models/common/swBoundaryConditionTimeSeries.h" // computeBC_E (time series)

// ---------------------------------------------------------------------------
//  Types AD (relogés depuis l'ancien swSp.h ; absents du natif).
// ---------------------------------------------------------------------------
struct SegPhysData {
    double   t_start;
    size_t   n_steps;
    std::vector<double> dt_seq;    // [step]
    std::vector<double> bc_flat;   // [step * edge_bnd_nb + j]
    std::vector<int>    obs_trig;  // [step]: -1 = no obs, >= 0 = obs index
    std::vector<double> h_ckpt, qx_ckpt, qy_ckpt;
};

// Changement de variable m -> phi appliqué DANS la tape Enzyme (input.solv 'phi_transform').
enum class PhiTransform {
    None = 0,   // phi = m
    Clamp,      // phi = clamp(m, 0, 1.05)
    Tanh,       // phi = 1/(1+exp(-lambda*m))
    Normalize   // phi = (m - lo)/(hi - lo)
};


class swSpE : public swSp
{
public:
    swSpE() : swSp() {}
    ~swSpE() = default;

    // =================================================================
    //  Helpers inline (relogés verbatim depuis l'ancien swSp.h)
    // =================================================================
    inline double phi_correct(double m) const {
        switch (phi_transform) {
            case PhiTransform::Clamp: {
                const double lo = 0.0, hi = 1.05;
                return (m < lo) ? lo : ((m > hi) ? hi : m);
            }
            case PhiTransform::Tanh:
                return 1.0 / (1.0 + std::exp(-lambda_tanh * m));
            case PhiTransform::Normalize:
                return (m - phi_norm_lo) * phi_norm_scale;
            case PhiTransform::None:
            default:
                return m;
        }
    }

    inline double dphi_dm(double m) const {
        switch (phi_transform) {
            case PhiTransform::Clamp: {
                const double lo = 0.0, hi = 1.05;
                return (m < lo || m > hi) ? 0.0 : 1.0;
            }
            case PhiTransform::Tanh: {
                const double s = 1.0 / (1.0 + std::exp(-lambda_tanh * m));
                return lambda_tanh * s * (1.0 - s);
            }
            case PhiTransform::Normalize:
                return phi_norm_scale;
            case PhiTransform::None:
            default:
                return 1.0;
        }
    }

    inline double m_from_phi(double phi) const {
        switch (phi_transform) {
            case PhiTransform::Tanh: {
                const double eps = 1e-2;
                const double p = (phi < eps) ? eps : ((phi > 1.0 - eps) ? 1.0 - eps : phi);
                return std::log(p / (1.0 - p)) / lambda_tanh;
            }
            case PhiTransform::Clamp:
                return (phi < 0.0) ? 0.0 : ((phi > 1.05) ? 1.05 : phi);
            case PhiTransform::Normalize:
            case PhiTransform::None:
            default:
                return phi;
        }
    }

    inline void set_quiet_logs(bool q) { minim_quiet_logs = q; }

    inline double get_adam_step_size()  const { return adam_step_size; }
    inline const std::string& get_minimize_algo() const { return minimize_algo; }
    inline double get_alpha_l1()  const { return alpha_l1; }
    inline double get_alpha_spatial() const { return alpha_spatial; }
    inline bool   get_zonation()  const { return zonation; }
    inline double get_alpha_phi_ref() const { return alpha_phi_ref; }
    inline double get_diag_grad_l2()    const { return diag_grad_l2; }
    inline double get_diag_grad_max()   const { return diag_grad_max; }
    inline long   get_diag_grad_argmax()const { return diag_grad_argmax; }
    inline double get_diag_dh_peak()    const { return diag_dh_peak; }

    // Accès porosité via phi_ptr (change de variable dans la tape).
    inline double& get_phi(size_t i) {
        return this->phi_ptr ? this->phi_ptr[i] : this->phi[i];
    }

    // Accesseur public au tableau natif `phi` (protected dans swSp) pour le
    // driver main_enzyme, qui n'est pas une méthode de swSpE.
    inline double* phi_data() { return phi.data(); }

    // =================================================================
    //  Boucle physique différenciée (ports Enzyme-safe)
    // =================================================================
    void   computeEdgeFactors_E();

    // interpolation temporelle Enzyme-safe (ex-swParent, inline). Membres BC de swParent
    // (bc_timeseries_list, bc_val_interpolated, edge_bnd_bc_pointer) tous accessibles.
    inline double interpolateTimeSeries_E(std::vector<double> const& times_vector,
                                          std::vector<double> const& vals_vector,
                                          double const& t)
    {
        size_t nb_times = times_vector.size();
        if (nb_times == 0) return 0.0;
        if (nb_times == 1 || t <= times_vector[0] - 1.0e-12) return vals_vector[0];
        if (t >= times_vector[nb_times - 1])                 return vals_vector[nb_times - 1];
        size_t iterator = 0;
        while (iterator < nb_times - 1 && t > times_vector[iterator + 1]) iterator++;
        double t_1 = times_vector[iterator],     t_2 = times_vector[iterator + 1];
        double v_1 = vals_vector[iterator],      v_2 = vals_vector[iterator + 1];
        double den = t_2 - t_1;
        if (den < 1e-12) return v_1;
        return v_1 + (v_2 - v_1) * (t - t_1) / den;
    }

    // computeBC_E : version SP (inline, ex-swSp.h) qui remplit bc_val_interpolated
    // INDEXÉ PAR ARÊTE (edge_bnd_nb). CRUCIAL : bc_type_nb vaut toujours 0 (jamais
    // initialisé), donc la version swParent (indexée par bc_type_nb) ne remplit RIEN
    // -> l'inflow n'est jamais injecté -> domaine sec. hypComputeFluxes_E lit
    // bc_val_interpolated[i] par arête i, d'où l'indexation par arête ici.
    inline void computeBC_E(double const& t) {
        for (size_t i = 0; i < edge_bnd_nb; i++) {
            swBoundaryConditionTimeSeries* bcts =
                static_cast<swBoundaryConditionTimeSeries*>(edge_bnd_bc_pointer[i]);
            bc_val_interpolated[i] = interpolateTimeSeries_E(
                bcts->time_values, bcts->bc_values, t);
        }
        for (swBoundaryConditionTimeSeries* bc : bc_timeseries_list) {
            bc->computeValBC(t);
        }
    }

    void   hypBoundaryFlux_E(const double& bR, const double& phiR, const double& zbR,
                             const double& zR, const double& qR, const double& rR,
                             const double& lambda_plus_R,
                             std::string const& bc_type, const double& bc_val,
                             double& flux1, double& flux2, double& flux3, double& source_right);
    void   hypCompute_E(size_t n_it_div_max, double& dt);
    void   hypComputeLRStates_E();
    double hypComputeTimeStep_E(void);
    void   hypComputeFluxes_E(void);
    void   hypComputeBalance_E(double dt);
    bool   hypMusclEvrComputeLRStates_E(void);
    bool   hypGodunovComputeLRStates_E(void);
    void   hypHllc_E(const double& bL, const double& bR,
                     const double& phiL, const double& phiR,
                     const double& zbL, const double& zbR,
                     const double& zL, const double& qL, const double& rL,
                     const double& zR, const double& qR, const double& rR,
                     const double& lambda_minus_L, const double& lambda_minus_R,
                     const double& lambda_plus_L, const double& lambda_plus_R,
                     double& flux1, double& flux2, double& flux3,
                     double& source_left, double& source_right);
    double specificForce_E(double boussinesq, double phi, double h, double u);
    void   sourceFrictionStrickler_E(double dt);
    void   sourceSingularHeadLoss_E(double dt);

    // =================================================================
    //  Pilotage AD / assimilation (checkpoint, gradient, minimisation)
    // =================================================================
    void   allocate_as_adjoint(swSpE* primal_solver);
    void   save_initial_state();
    void   reset_to_start();
    double compute_cost_function();
    double run_pure_simulation(double K_val);
    double run_pure_simulation(double* phi_map);
    void   test_enzyme_integration(swSpE* solver);
    void   compare_versions(double phi_val);
    void   init_primal_memory();
    double compute_cost_frame_ptr(const double* h_obs);
    double compute_cost_frame(const std::vector<double>& h_obs);

    bool   load_all_observations();
    bool   load_single_observation_file(double time, const std::string& filepath);
    bool   load_phi_csv(const std::string& filepath, double* phi_out);

    double execute_E(double* phi_map);
    double execute_E_for_enzyme(double* phi_map, const double* bc_precomp, int n_bc_steps);

    void   get_spatial_gradient_porosity(double* phi_map, double* dJ_dPhi);
    void   get_spatial_gradient_checkpointed(double* phi_map, double* dJ_dPhi,
                                             size_t checkpoint_stride = 30);
    double replay_trajectory_cost(double* phi_map);
    void   run_taylor_test_checkpointed(double* phi_base, double* dJ_dPhi,
                                        size_t n_top_cells = 3);
    void   run_taylor_test_porosity(double* phi_base, double* gradient_enzyme);

    bool   begin_reference_generation();
    bool   finalize_reference_generation();
    bool   generate_references_E(double* phi_map, size_t checkpoint_stride = 30);
    bool   run_forward_E(double* phi_map, size_t checkpoint_stride = 30);

    double run_minimization(double* phi_start, size_t max_it, size_t checkpoint_stride,
                            double alpha_reg);
    void   build_checkpoint_records(double* phi_map, size_t stride);

    // Fonction statique différenciée par Enzyme (prend swSpE* — cf. règle C++ ci-dessus).
    static double ckpt_segment_E(
        swSpE*       solver,
        double*      phi_map,
        const double* h_in,
        const double* qx_in,
        const double* qy_in,
        const double* dt_seq,
        const double* bc_flat,
        const int*    obs_trig,
        size_t        n_steps,
        size_t        edge_bnd_nb,
        const double* obs_times,
        const double* const* obs_frames,
        size_t        n_obs,
        double        t_start);

    void   update_phi_norm_bounds(const double* m_map);
    void   resolve_cost_weight(const std::string& spec, std::vector<double>& out, const char* name);
    void   storeReferenceMap(const double& t);   // AD : dump carte d'état (gen-ref / forward dump)

    // ---- Assimilation aux sondes (obs_type=="gauges"), chemin séparé ----
    bool   locate_gauges(swMesh* mesh);
    bool   load_gauge_observations(const std::string& file);
    double compute_cost_gauges_ptr(const double* obs_g);

    // =================================================================
    //  Membres de données AD (classés : présents dans le courant, absents
    //  du natif). Les membres natifs (phi, water_depth, ...) sont hérités.
    // =================================================================
    // -- porosité / change de variable
    double* phi_ptr = nullptr;
    PhiTransform phi_transform = PhiTransform::None;
    double  lambda_tanh   = 1.0;
    double  phi_norm_lo   = 0.0;
    double  phi_norm_hi   = 1.0;
    double  phi_norm_scale = 1.0;
    std::vector<size_t> phi_norm_free_idx;
    double  phi_lower_bound = 0.0;
    double  phi_upper_bound = 1.0;

    // -- buffers reconstruction / hyperbolique (membres, ex-locaux natifs)
    std::vector<double> z_cell = {};
    std::vector<double> dz_dx  = {};
    std::vector<double> dz_dy  = {};
    std::vector<double> z_slope_x = {}, z_slope_y = {};
    std::vector<double> q_slope_x = {}, q_slope_y = {};
    std::vector<double> r_slope_x = {}, r_slope_y = {};
    std::vector<double> swept_data;

    // -- état initial sauvegardé (reset_to_start)
    std::vector<double> h_init = {}, qx_init = {}, qy_init = {};

    // -- checkpoints / trajectoire préenregistrée
    std::vector<SegPhysData> segs_record_;
    std::vector<double> dt_seq;    // [step]
    std::vector<double> bc_flat;   // [step * edge_bnd_nb + j]
    std::vector<int>    obs_trig;  // [step]
    std::vector<double> h_ckpt, qx_ckpt, qy_ckpt;
    int    adjoint_truncation_seg_ = -1;
    int    forced_first_seg_ = -1;
    std::vector<size_t> assim_window_tiles_;

    // -- observations
    std::map<double, std::vector<double>> observations_map;
    std::vector<double> obs_times;
    std::vector<double> obs_times_vec;
    std::vector<std::vector<double>> obs_value;
    std::vector<std::vector<double>> obs_frames_storage;
    std::vector<const double*> obs_frames_ptrs;
    const double* current_obs_times = nullptr;
    const double* const* current_obs_frames = nullptr;
    size_t current_n_obs = 0;
    std::vector<double> water_depth_obs = {};

    // -- génération de références / dump
    bool   gen_ref_mode = false;
    size_t n_ref_written = 0;
    std::string state_dump_dir = "reference";
    bool   forward_dump_mode = false;

    // -- régularisation / optimiseur
    std::string minimize_algo = "bfgs";
    double adam_step_size = 0.05;
    bool   adam_safe = false;
    double alpha_l1 = 0.0;
    double alpha_phi_ref = 0.0;
    double alpha_spatial = 0.0;
    int    lbfgsb_m = 6;
    double lbfgsb_epsilon = 1e-5;
    int    lbfgsb_max_linesearch = 20;

    // -- fenêtre d'assimilation / pavage
    double t_state_assimilation = 0.0;
    double dt_map_ad = -1.0;
    double assimilation_window_duration = -1.0;
    bool   windowed_assimilation = false;
    long   assimilation_window_steps = -1;
    long   n_assimilation_windows = 1;
    bool   zonation = false;
    bool   mask_dry_cells = false;
    bool   mask_wall_gradient = false;
    bool   cost_normalize = false;

    // -- pondération du coût
    std::string h_weight_spec = "1.0";
    std::string q_weight_spec = "1.0";
    std::vector<double> h_weight_map;   // taille cell_nb
    std::vector<double> q_weight_map;   // taille cell_nb

    // -- chemins
    std::string porosity_ref_map_path;
    std::string porosity_free_zone_map_path;

    // -- diagnostics de descente
    double diag_grad_l2     = 0.0;
    double diag_grad_max    = 0.0;
    long   diag_grad_argmax = -1;
    double diag_dh_peak     = 0.0;
    bool   taylor_test_after_gradient = true;
    bool   minim_quiet_logs = false;

    // -- assimilation aux sondes (obs_type=="gauges")
    std::string obs_type = "cells";
    size_t n_gauges = 0;
    std::vector<size_t> gauge_cell_idx;
    std::vector<size_t> gauge_src_col;
    std::vector<double> gauge_w;
};

// Macro AD : accès porosité via phi_ptr. N'expanse QUE la syntaxe d'appel phi(i)
// (jamais phi[i] ni la déclaration `std::vector<double> phi`). Confinée aux unités
// de compilation assim/. Définie APRÈS l'inclusion du natif pour ne pas le perturber.
#define phi(i) (this->phi_ptr ? this->phi_ptr[i] : this->phi[i])
