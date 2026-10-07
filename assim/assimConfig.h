#pragma once
// =====================================================================
//  assimConfig — lecture des clés d'assimilation/AD depuis input.solv,
//  DÉCOUPLÉE du parseur natif swSolverParam.
//
//  Le parseur natif pristine `logFatal` sur toute clé inconnue. On lit donc
//  ICI les ~32 clés AD, et on écrit une COPIE FILTRÉE d'input.solv (clés AD
//  retirées) dans un dossier temporaire que le parseur natif
//  swSolverParam::readParameters() lit tel quel. Zéro modif native ;
//  réutilisation maximale de la fonction de lecture native.
// =====================================================================
#include <string>
#include <cstddef>

class swSpE;

struct AssimConfig {
    // -- Clés poussées dans les membres de swSpE (apply_to) --
    std::string minimize_algo = "bfgs";
    double adam_step_size = 0.05;
    bool   adam_safe = false;
    double alpha_l1 = 0.0;
    double alpha_phi_ref = 0.0;
    double alpha_spatial = 0.0;
    double lambda_tanh = 1.0;
    double lbfgsb_epsilon = 1e-5;
    int    lbfgsb_m = 6;
    int    lbfgsb_max_linesearch = 20;
    std::string phi_transform = "none";      // none|clamp|tanh|normalize
    double phi_lower_bound = 0.0;
    double phi_upper_bound = 1.0;
    double t_state_assimilation = 0.0;
    double dt_map_ad = -1.0;
    bool   windowed_assimilation = false;
    double assimilation_window_duration = -1.0;
    long   assimilation_window_steps = -1;
    long   n_assimilation_windows = 1;
    bool   zonation = false;
    bool   mask_dry_cells = false;
    bool   mask_wall_gradient = false;
    bool   cost_normalize = false;
    std::string obs_type = "cells";
    std::string h_weight = "1.0";
    std::string q_weight = "1.0";
    std::string porosity_ref_map;
    std::string porosity_free_zone_map;

    // -- Clés utilisées directement par main_enzyme --
    double alpha_reg = 0.0;
    std::size_t checkpoint_stride = 30;
    std::size_t max_it = 3;
    std::string porosity_ebauche_map;

    // Lit les clés AD de `in_solv_path` (input.solv utilisateur) et écrit dans
    // `filtered_out_path` une copie SANS les clés AD (commentaires inline retirés
    // sur les lignes de clés natives), destinée au parseur natif. Retourne false
    // si un fichier ne peut être ouvert.
    bool load_and_filter(const std::string& in_solv_path,
                         const std::string& filtered_out_path);

    // Pousse les valeurs AD dans les membres publics de swSpE (après le
    // importAttributesFromSolverParam natif, qui ne connaît pas ces clés).
    void apply_to(swSpE* m) const;
};
