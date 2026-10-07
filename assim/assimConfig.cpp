// =====================================================================
//  assimConfig.cpp — implémentation du lecteur de config AD découplé.
// =====================================================================
#include "assimConfig.h"
#include "swSpE.h"

#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <cctype>

namespace {

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

bool parse_bool(const std::string& v) {
    std::string l = to_lower(v);
    return (l == "1" || l == "true" || l == "yes" || l == "on");
}

// Ensemble figé des clés AD (ignorées par le parseur natif -> filtrées ici).
const std::set<std::string>& ad_keys() {
    static const std::set<std::string> K = {
        "adam_safe", "adam_step_size", "alpha_l1", "alpha_phi_ref", "alpha_reg",
        "alpha_spatial", "assimilation_window_duration", "assimilation_window_steps",
        "checkpoint_stride", "cost_normalize", "dt_map_ad", "h_weight", "lambda_tanh",
        "lbfgsb_epsilon", "lbfgsb_m", "lbfgsb_max_linesearch", "mask_dry_cells",
        "mask_wall_gradient", "max_it", "minimize_algo", "n_assimilation_windows",
        "obs_type", "phi_lower_bound", "phi_transform", "phi_upper_bound",
        "porosity_ebauche_map", "porosity_free_zone_map", "porosity_ref_map",
        "q_weight", "t_state_assimilation", "windowed_assimilation", "zonation"
    };
    return K;
}

} // namespace

void AssimConfig::apply_to(swSpE* m) const {
    m->minimize_algo        = minimize_algo;
    m->adam_step_size       = adam_step_size;
    m->adam_safe            = adam_safe;
    m->alpha_l1             = alpha_l1;
    m->alpha_phi_ref        = alpha_phi_ref;
    m->alpha_spatial        = alpha_spatial;
    m->lambda_tanh          = lambda_tanh;
    m->lbfgsb_epsilon       = lbfgsb_epsilon;
    m->lbfgsb_m             = lbfgsb_m;
    m->lbfgsb_max_linesearch= lbfgsb_max_linesearch;
    m->phi_lower_bound      = phi_lower_bound;
    m->phi_upper_bound      = phi_upper_bound;
    m->t_state_assimilation = t_state_assimilation;
    m->dt_map_ad            = dt_map_ad;
    m->windowed_assimilation= windowed_assimilation;
    m->assimilation_window_duration = assimilation_window_duration;
    m->assimilation_window_steps    = assimilation_window_steps;
    m->n_assimilation_windows       = n_assimilation_windows;
    m->zonation             = zonation;
    m->mask_dry_cells       = mask_dry_cells;
    m->mask_wall_gradient   = mask_wall_gradient;
    m->cost_normalize       = cost_normalize;
    m->obs_type             = obs_type;
    m->h_weight_spec        = h_weight;
    m->q_weight_spec        = q_weight;
    m->porosity_ref_map_path       = porosity_ref_map;
    m->porosity_free_zone_map_path = porosity_free_zone_map;

    // string -> enum PhiTransform (défini dans swSpE.h)
    const std::string pt = to_lower(phi_transform);
    if      (pt == "clamp")     m->phi_transform = PhiTransform::Clamp;
    else if (pt == "tanh")      m->phi_transform = PhiTransform::Tanh;
    else if (pt == "normalize") m->phi_transform = PhiTransform::Normalize;
    else                        m->phi_transform = PhiTransform::None;
}

bool AssimConfig::load_and_filter(const std::string& in_solv_path,
                                  const std::string& filtered_out_path) {
    std::ifstream fin(in_solv_path);
    if (!fin) return false;
    std::ofstream fout(filtered_out_path);
    if (!fout) return false;

    const auto& AD = ad_keys();
    std::string line;
    while (std::getline(fin, line)) {
        // Lignes commentaire/blanches : le natif les saute (premier char # ou =).
        // On les recopie telles quelles.
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) { fout << line << "\n"; continue; }
        char c0 = line[first];
        if (c0 == '#' || c0 == '=') { fout << line << "\n"; continue; }

        // Ligne de contenu : retirer un éventuel commentaire inline (# ...).
        std::string content = line;
        size_t hash = content.find('#');
        if (hash != std::string::npos) content = content.substr(0, hash);

        std::istringstream iss(content);
        std::string key, val;
        iss >> key >> val;
        if (key.empty()) { continue; }                 // ligne vide après nettoyage
        std::string lkey = to_lower(key);

        if (AD.count(lkey)) {
            // Clé AD : capturée ici, PAS écrite dans la copie native.
            if      (lkey == "minimize_algo")        minimize_algo = val;
            else if (lkey == "adam_step_size")       adam_step_size = std::stod(val);
            else if (lkey == "adam_safe")            adam_safe = parse_bool(val);
            else if (lkey == "alpha_l1")             alpha_l1 = std::stod(val);
            else if (lkey == "alpha_phi_ref")        alpha_phi_ref = std::stod(val);
            else if (lkey == "alpha_reg")            alpha_reg = std::stod(val);
            else if (lkey == "alpha_spatial")        alpha_spatial = std::stod(val);
            else if (lkey == "lambda_tanh")          lambda_tanh = std::stod(val);
            else if (lkey == "lbfgsb_epsilon")       lbfgsb_epsilon = std::stod(val);
            else if (lkey == "lbfgsb_m")             lbfgsb_m = std::stoi(val);
            else if (lkey == "lbfgsb_max_linesearch")lbfgsb_max_linesearch = std::stoi(val);
            else if (lkey == "phi_transform")        phi_transform = val;
            else if (lkey == "phi_lower_bound")      phi_lower_bound = std::stod(val);
            else if (lkey == "phi_upper_bound")      phi_upper_bound = std::stod(val);
            else if (lkey == "t_state_assimilation") t_state_assimilation = std::stod(val);
            else if (lkey == "dt_map_ad")            dt_map_ad = std::stod(val);
            else if (lkey == "windowed_assimilation")windowed_assimilation = parse_bool(val);
            else if (lkey == "assimilation_window_duration") assimilation_window_duration = std::stod(val);
            else if (lkey == "assimilation_window_steps")    assimilation_window_steps = std::stol(val);
            else if (lkey == "n_assimilation_windows")       n_assimilation_windows = std::stol(val);
            else if (lkey == "zonation")             zonation = parse_bool(val);
            else if (lkey == "mask_dry_cells")       mask_dry_cells = parse_bool(val);
            else if (lkey == "mask_wall_gradient")   mask_wall_gradient = parse_bool(val);
            else if (lkey == "cost_normalize")       cost_normalize = parse_bool(val);
            else if (lkey == "obs_type")             obs_type = val;
            else if (lkey == "h_weight")             h_weight = val;
            else if (lkey == "q_weight")             q_weight = val;
            else if (lkey == "porosity_ref_map")     porosity_ref_map = val;
            else if (lkey == "porosity_free_zone_map") porosity_free_zone_map = val;
            else if (lkey == "porosity_ebauche_map") porosity_ebauche_map = val;
            else if (lkey == "checkpoint_stride")    checkpoint_stride = (std::size_t)std::stoul(val);
            else if (lkey == "max_it")               max_it = (std::size_t)std::stoul(val);
        } else {
            // Clé native : écrite proprement (clé + valeur, sans commentaire inline)
            // -> le parseur natif exige exactement 2 tokens.
            fout << key << " " << val << "\n";
        }
    }
    return true;
}
