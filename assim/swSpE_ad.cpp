#include "swSpE.h"
#include "../swParent.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>
#include "../../../common/swLog.h"
#include "../../../common/swMessage.h"
#include "../../../common/swIO.h"
#include <cmath>

// --- Constantes Enzyme ---
extern int enzyme_dup;
extern int enzyme_const;
extern int enzyme_out;
extern int enzyme_dupnoneed; 

    


/**
 * @brief Fonction de coût (Objectif)
 * Calcule l'écart quadratique entre la hauteur d'eau simulée et observée.
 */


/*double swSpE::compute_cost_frame(const std::vector<double>& h_obs) { // fonction cout
    
    //return 47.0;
    double J = 0.0;
    // Calcul de l'écart quadratique sur les cellules
    for (size_t i = 0; i < cell_nb; i++) {
        double diff = water_depth[i] - h_obs[i];
        J += diff * diff;
    }
    // On normalise par le nombre de cellules
    return J / static_cast<double>(cell_nb);
}*/




bool swSpE::load_all_observations()
{
    namespace fs = std::filesystem;

    // Dossier contenant les observations de reference : <dossier_simu>/reference/
    const std::string ref_dir = m_io->absolutePath("reference");

    if (!fs::exists(ref_dir) || !fs::is_directory(ref_dir)) {
        logError() << "Dossier d'observations introuvable : " << ref_dir;
        logError() << "  Les observations doivent etre des fichiers reference/ResultHydro_p*.txt";
        return false;
    }

    // Fenetre d'assimilation : [t_state_assimilation, t_max], bornes lues depuis input.solv
    const double t_begin = t_state_assimilation;
    const double t_end   = t_max;
    const double eps     = 1e-6; // tolerance sur les bornes de la fenetre

    logInfo() << "Chargement des observations dans la fenetre d'assimilation ["
              << t_begin << "s, " << t_end << "s] depuis " << ref_dir;

    observations_map.clear();

    // On collecte d'abord les fichiers candidats pour pouvoir les trier/journaliser
    std::vector<std::pair<double, std::string>> candidates; // (temps, chemin_complet)
    int scanned = 0;
    int skipped = 0;

    for (const auto& entry : fs::directory_iterator(ref_dir)) {
        if (!entry.is_regular_file()) continue;

        const std::string fname = entry.path().filename().string();

        // On ne garde que les fichiers du type ResultHydro_p*.txt
        if (fname.rfind("ResultHydro_p", 0) != 0) continue;       // ne commence pas par le prefixe
        if (entry.path().extension() != ".txt") continue;

        scanned++;

        // Le temps (en secondes) est encode dans le nom du fichier
        const double t_obs = swIO::timeFromFileName(fname);

        // Filtrage sur la fenetre d'assimilation
        if (t_obs < t_begin - eps || t_obs > t_end + eps) {
            skipped++;
            continue;
        }

        candidates.emplace_back(t_obs, entry.path().string());
    }

    // Tri par temps croissant (la map est de toute facon ordonnee, mais on journalise dans l'ordre)
    std::sort(candidates.begin(), candidates.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    int loaded = 0;
    for (const auto& [t_obs, full_path] : candidates) {
        std::cout << "[OBS] t=" << t_obs << "s  ->  " << full_path << std::endl;

        if (!load_single_observation_file(t_obs, full_path)) {
            logError() << "Echec du chargement de l'observation t=" << t_obs
                       << "  fichier : " << full_path;
            return false;
        }
        loaded++;
    }

    if (loaded == 0) {
        logError() << "Aucune observation ResultHydro_p*.txt trouvee dans " << ref_dir
                   << " pour la fenetre [" << t_begin << "s, " << t_end << "s]";
        logError() << "  (" << scanned << " fichier(s) scanne(s), " << skipped
                   << " hors fenetre)";
        return false;
    }

    logInfo() << loaded << " observation(s) chargee(s) depuis " << ref_dir
              << " (" << skipped << " hors fenetre, " << scanned << " scanne(s))";
    return true;
}

bool swSpE::load_single_observation_file(double time, const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) return false;

    // On prépare un vecteur pour stocker h, qx et qy (taille x 3)
    std::vector<double> current_frame(cell_nb * 3, 0.0);
    std::string line;
    std::getline(file, line); // Sauter l'en-tête

    int cell_count = 0;
    while (std::getline(file, line) && cell_count < cell_nb) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        double val;
        
        // On lit les colonnes : 0:x, 1:y, 2:zb, 3:h, 4:qx, 5:qy
        for (int col = 0; col <= 5; ++col) {
            if (!(ss >> val)) break;
            if (col == 3) current_frame[cell_count] = val;               // h
            if (col == 4) current_frame[cell_count + cell_nb] = val;     // qx
            if (col == 5) current_frame[cell_count + 2 * cell_nb] = val; // qy
        }
        cell_count++;
    }
    file.close();
    observations_map[time] = current_frame;
    return true;
}

void swSpE::init_primal_memory() {
    // Redimensionnement des tableaux de TRAVAIL uniquement
    // NE PAS toucher à water_depth, unit_discharge_x/y — ce sont les CI
    
    this->z_cell.assign(this->cell_nb, 0.0);
    this->dz_dx.assign(this->cell_nb, 0.0);
    this->dz_dy.assign(this->cell_nb, 0.0);

    this->z_slope_x.assign(this->cell_nb, 0.0);
    this->z_slope_y.assign(this->cell_nb, 0.0);
    this->q_slope_x.assign(this->cell_nb, 0.0);
    this->q_slope_y.assign(this->cell_nb, 0.0);
    this->r_slope_x.assign(this->cell_nb, 0.0);
    this->r_slope_y.assign(this->cell_nb, 0.0);

    this->swept_data.assign(this->cell_nb, 0.0);

    // Friction — redimensionnement si nécessaire, valeurs par défaut
    if (this->friction_coeff_direction_1.size() != this->cell_nb)
        this->friction_coeff_direction_1.assign(this->cell_nb, 0.0);
    if (this->friction_coeff_direction_2.size() != this->cell_nb)
        this->friction_coeff_direction_2.assign(this->cell_nb, 0.0);
    if (this->friction_cos_alpha.size() != this->cell_nb)
        this->friction_cos_alpha.assign(this->cell_nb, 1.0);
    if (this->friction_sin_alpha.size() != this->cell_nb)
        this->friction_sin_alpha.assign(this->cell_nb, 0.0);

        
    // Interfaces internes
    this->h_riemann_int.assign(edge_int_nb, {0.0, 0.0});
    this->z_riemann_int.assign(edge_int_nb, {0.0, 0.0});
    this->q_riemann_int.assign(edge_int_nb, {0.0, 0.0});
    this->r_riemann_int.assign(edge_int_nb, {0.0, 0.0});
    this->lambda_minus_int.assign(edge_int_nb, {0.0, 0.0});
    this->lambda_plus_int.assign(edge_int_nb, {0.0, 0.0});
    this->lambda_max_int.assign(edge_int_nb, 0.0);

    // Interfaces de bord
    this->h_riemann_bnd.assign(edge_bnd_nb, 0.0);
    this->z_riemann_bnd.assign(edge_bnd_nb, 0.0);
    this->q_riemann_bnd.assign(edge_bnd_nb, 0.0);
    this->r_riemann_bnd.assign(edge_bnd_nb, 0.0);
    this->lambda_minus_bnd.assign(edge_bnd_nb, 0.0);
    this->lambda_plus_bnd.assign(edge_bnd_nb, 0.0);
    this->lambda_max_bnd.assign(edge_bnd_nb, 0.0);

    // phi : redimensionner seulement si nécessaire, NE PAS écraser
    if (this->phi.size() != this->cell_nb)
        this->phi.assign(this->cell_nb, 1.0);

    // water_depth et unit_discharge : NE PAS TOUCHER
    // Elles contiennent les CI réelles chargées par setInitialConditions()
    if (this->water_depth.size() != this->cell_nb)
        this->water_depth.assign(this->cell_nb, 0.0);
    if (this->unit_discharge_x.size() != this->cell_nb)
        this->unit_discharge_x.assign(this->cell_nb, 0.0);
    if (this->unit_discharge_y.size() != this->cell_nb)
        this->unit_discharge_y.assign(this->cell_nb, 0.0);
}

double swSpE::execute_E(double* phi_map) {
    if (this->cell_nb <= 0) return 1e30;

    this->reset_to_start();
    
    this->phi_ptr = phi_map;
    for (size_t i = 0; i < this->cell_nb; i++) {
        this->phi[i] = phi_map[i];
    }

    double total_accumulated_cost = 0.0;
    size_t obs_idx = 0; 

    while (t < t_max) {
        computeBC_E(t);

        double dt_target = dt_max;
        bool trigger_obs = false;
        
        if (obs_idx < this->current_n_obs) {
            double t_next = this->current_obs_times[obs_idx];
            dt_target = t_next - t; 
        }

        double current_dt = dt_target;
        
        // hypCompute_E va maintenant modifier proprement current_dt si la CFL est plus restrictive
        hypCompute_E(n_it_div_max, current_dt);

        // Protection contre les miettes d'arrondis machine
        if (obs_idx < this->current_n_obs && std::abs(current_dt - dt_target) < 1e-11) {
            current_dt = dt_target; 
            trigger_obs = true;
        }

        sourceFrictionStrickler_E(current_dt);
        // Head loss singulier (gardé) : cohérence avec build_checkpoint_records /
        // ckpt_segment_E, sinon J_execute (contrôle du test de Taylor) diverge de J_traj.
        if (model_singular_head_loss != "none")
            sourceSingularHeadLoss_E(current_dt);

        if (std::isnan(this->water_depth[0])) { return 1e30; }

        t += current_dt;

        if (trigger_obs) {
            total_accumulated_cost += (obs_type == "gauges")
                ? this->compute_cost_gauges_ptr(this->current_obs_frames[obs_idx])
                : this->compute_cost_frame_ptr(this->current_obs_frames[obs_idx]);
            obs_idx++;
        }
    }

    return total_accumulated_cost;
}


double swSpE::compute_cost_frame_ptr(const double* obs_ptr) {
    double J = 0.0;
    
    // Décalages pour accéder aux différentes variables dans le bloc
    const double* h_obs  = obs_ptr;
    const double* qx_obs = obs_ptr + cell_nb;
    const double* qy_obs = obs_ptr + 2 * cell_nb;

    // Ponderation spatiale par cellule (h_weight_map / q_weight_map, resolues dans init()).
    // Cartes vides (securite) -> poids 1. NB : ces cartes ne derivent PAS de phi -> Enzyme
    // les voit comme des constantes (la fonction est differenciee dans ckpt_segment_E).
    const bool wh_ok = (h_weight_map.size() == cell_nb);
    const bool wq_ok = (q_weight_map.size() == cell_nb);
    for (size_t i = 0; i < cell_nb; i++) {
        // Écart sur la hauteur d'eau
        double dh = water_depth[i] - h_obs[i];
        // Écart sur les débits unitaires
        double dqx = unit_discharge_x[i] - qx_obs[i];
        double dqy = unit_discharge_y[i] - qy_obs[i];

        const double wh = wh_ok ? h_weight_map[i] : 1.0;
        const double wq = wq_ok ? q_weight_map[i] : 1.0;

        // Coût L2 pondéré : wh*dh^2 + wq*(dqx^2 + dqy^2)
        J += wh * (dh * dh) + wq * ((dqx * dqx) + (dqy * dqy));
    }
    // Normalisation optionnelle (input.solv 'cost_normalize') : coût MOYEN par
    // observation au lieu de la SOMME -> echelle de J (et du gradient) invariante a
    // la largeur de la fenetre d'assimilation. current_n_obs ne derive pas de phi
    // (constante pour Enzyme) -> le gradient est divise par le meme facteur, donc
    // valeur et gradient restent coherents (et le test de Taylor reste a 1).
    if (cost_normalize && current_n_obs > 0)
        J /= (double)current_n_obs;
    return J;
    //return this->compute_cost_frame_ptr(this->current_obs_frames[this->current_n_obs - 1]);
}

double swSpE::compute_cost_frame(const std::vector<double>& h_obs) {
    return compute_cost_frame_ptr(h_obs.data());
}


// --- Wrappers pour Enzyme ---

struct EnzymeObsContext {
    const double* times;
    const double* const* frames;
    size_t n_obs;
};

//=== TEST Changment appel a enzyme

/*static void enzyme_wrapper_porosity(swSpE* solver, double* phi_map, double* J_out) {
    // deuxième argument double* phi_map, suprimé , apéle a enzuyme modifier , on vas voire comment enzyme s'en sort.
    // Le solver contient déjà les pointeurs vers les observations
    *J_out = solver->execute_E(phi_map);
}*/

// [UNIFICATION ENZYME] Wrapper et extern __enzyme_autodiff NEUTRALISES.
// Ce wrapper differenciait execute_E (sans computeEdgeFactors_E -> bug C).
// Plus aucun __enzyme_autodiff dans ce fichier : le seul execute differencie
// du projet est ckpt_segment_E (porosity_checkpoint.cpp). Conserve en
// commentaire pour memoire historique.
//
// double enzyme_wrapper_porosity(swSpE* solver, double* phi_array) {
//     return solver->execute_E(phi_array);
// }
//
// extern "C" {
//     void __enzyme_autodiff(void*, ...);
// }


void swSpE::get_spatial_gradient_porosity(double* phi_map, double* dJ_dPhi)
{
    // =====================================================================
    // [UNIFICATION ENZYME] Chemin historique NEUTRALISE.
    //
    // Avant : cette fonction differenciait execute_E via __enzyme_autodiff.
    // Or hypCompute_E (appele par execute_E) n'appelle PAS computeEdgeFactors_E
    // -> le terme 1/phi restait HORS de la tape (bug C) -> gradient faux/nul
    // sur les cellules-barrage.
    //
    // Desormais on delegue au SEUL execute differencie du projet :
    //   ckpt_segment_E (computeEdgeFactors_E DANS la tape), pilote par
    //   get_spatial_gradient_checkpointed.
    // => un unique __enzyme_autodiff pour la porosite, gradient coherent,
    //    et coherent avec le forward replay_trajectory_cost (meme fonction).
    // =====================================================================

    const size_t n_obs = observations_map.size();
    if (n_obs == 0) {
        logError() << "get_spatial_gradient_porosity : aucune observation chargee "
                      "(observations_map vide) -> gradient nul.";
        std::fill(dJ_dPhi, dJ_dPhi + this->cell_nb, 0.0);
        return;
    }

    // Cablage des observations attendu par le chemin checkpointe
    // (memes etapes que main_enzyme avant get_spatial_gradient_checkpointed).
    this->obs_times_vec.clear();
    this->obs_frames_storage.clear();
    this->obs_frames_ptrs.clear();
    this->obs_frames_storage.reserve(n_obs);
    for (auto const& [t_obs, frame] : observations_map) {
        this->obs_times_vec.push_back(t_obs);
        this->obs_frames_storage.push_back(frame);
    }
    for (size_t i = 0; i < n_obs; ++i)
        this->obs_frames_ptrs.push_back(this->obs_frames_storage[i].data());
    this->current_obs_times  = this->obs_times_vec.data();
    this->current_obs_frames = this->obs_frames_ptrs.data();
    this->current_n_obs      = n_obs;

    // Delegation : stride par defaut (cf. signature checkpointed = 30).
    this->get_spatial_gradient_checkpointed(phi_map, dJ_dPhi);
}


void swSpE::run_taylor_test_porosity(double* phi_base, double* gradient_enzyme)
{
    std::cout << "=== DÉBUT TAYLOR DIAGNOSTIC ===" << std::endl;
    const size_t N = this->cell_nb;

    // 1. Direction de perturbation uniforme
    std::vector<double> delta_phi(N, 0.01);

    // 2. Produit scalaire : <∇J, δφ>
    double dot_product = 0.0;
    for (size_t i = 0; i < N; i++) {
        dot_product += gradient_enzyme[i] * delta_phi[i];
    }

    // ============================================================
    // CORRECTION CRITIQUE : Warm-up pour éliminer l'hystérésis d'état
    // ============================================================
    this->reset_to_start();
    this->execute_E(phi_base); 

    // Calcul du VRAI J_base dans cet état stabilisé
    this->reset_to_start();
    double J_base = this->execute_E(phi_base);
    // ============================================================

    std::cout << "J(base) stabilisé                     = " << J_base << std::endl;
    std::cout << "Produit scalaire <Grad_Enzyme, d_phi> = " << dot_product << std::endl;
    
    // --- OUVERTURE DU FICHIER CSV POUR L'EXPORT ---
    std::ofstream csv_file("taylor_test_porosity.csv");
    csv_file << "epsilon,J_eps,ratio,abs_error,gradient_fd,dot_product\n";
    csv_file << std::scientific << std::setprecision(10);
    // -----------------------------------------------

    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << " epsilon    | J(phi+eps*dphi) | Ratio Taylor | Erreur Abs." << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

    double epsilon = 0.10;
    std::vector<double> phi_perturbed(N);

    // 3. Boucle de convergence
    for (int k = 0; k < 7; k++) {
        for (size_t i = 0; i < N; i++) {
            phi_perturbed[i] = phi_base[i] + epsilon * delta_phi[i];
        }
        
        this->reset_to_start();
        double J_eps = this->execute_E(phi_perturbed.data());

        // Calcul du gradient local propre
        double gradient_fd = (J_eps - J_base) / epsilon;
        
        double ratio = (std::abs(dot_product) > 1e-15) ? (gradient_fd / dot_product) : 0.0;
        double error = std::abs(gradient_fd - dot_product);

        // Affichage console
        printf("%10.1e | %15.8e | %12.6f | %12.6e\n", epsilon, J_eps, ratio, error);

        // Écriture dans le CSV
        csv_file << epsilon << "," 
                 << J_eps << "," 
                 << ratio << "," 
                 << error << "," 
                 << gradient_fd << "," 
                 << dot_product << "\n";

        epsilon *= 0.1;
    }
    csv_file.close();
    std::cout << "============================================================\n" << std::endl;
    std::cout << "[Info] Données Taylor exportées dans 'taylor_test_porosity.csv'" << std::endl;

    // 4. Gradient FD par cellule (bénéficie aussi de l'état stabilisé)
    const double eps_fd = 1e-5;
    std::cout << "\n=== GRADIENT FD PAR CELLULE (eps=" << eps_fd << ") ===" << std::endl;
    
    this->reset_to_start();
    std::vector<double> phi_base_copy(phi_base, phi_base + N);
    double J0 = this->execute_E(phi_base_copy.data());

    for (size_t i = 0; i < N; i++) {
        std::vector<double> phi_p(phi_base, phi_base + N);
        phi_p[i] += eps_fd;
        
        this->reset_to_start();
        double Jp = this->execute_E(phi_p.data());
        double grad_fd_i = (Jp - J0) / eps_fd;
        
        printf("  dJ/dphi[%zu]  FD = %12.6f  Enzyme = %12.6f  ratio = %12.6f\n",
               i, grad_fd_i, gradient_enzyme[i],
               (std::abs(grad_fd_i) > 1e-15) ? gradient_enzyme[i] / grad_fd_i : 0.0);
    }
    std::cout << "=== FIN GRADIENT FD ===" << std::endl;
}

void swSpE::compare_versions(double phi_val) {
    // --- 1. SÉCURISATION ET INITIALISATION DES OBSERVATIONS ---
    // Obligatoire pour éviter que execute_E ne lise des pointeurs non initialisés (nullptr)
    size_t n_obs = this->observations_map.size();
    if (n_obs == 0) {
        std::cout << "[ERROR Compare] Aucune observation n'est chargée dans la map !" << std::endl;
        return;
    }

    this->obs_times_vec.clear();
    this->obs_frames_storage.clear();
    this->obs_frames_ptrs.clear();
    this->obs_frames_storage.reserve(n_obs);

    for (auto const& [t_obs, frame] : this->observations_map) {
        this->obs_times_vec.push_back(t_obs);
        this->obs_frames_storage.push_back(frame);
    }
    for (size_t i = 0; i < n_obs; i++) {
        this->obs_frames_ptrs.push_back(this->obs_frames_storage[i].data());
    }

    // On ancre solidement les pointeurs membres de la classe pour la durée du test
    this->current_obs_times  = this->obs_times_vec.data();
    this->current_obs_frames = this->obs_frames_ptrs.data();
    this->current_n_obs      = n_obs;


    // --- 2. PRÉPARATION DU TABLEAU DE TEST ---
    std::vector<double> test_phi(this->cell_nb, phi_val);
    double* old_phi_ptr = this->phi_ptr; // Sauvegarde préventive


    // --- 3. TEST VERSION ENZYME ---
    this->reset_to_start();
    this->phi_ptr = test_phi.data(); // On oriente le pointeur brut vers nos données de test
    
    double J_E = this->execute_E(test_phi.data());


    // --- 4. TEST VERSION ORIGINALE ---
    this->reset_to_start();
    
    // Crucial : On coupe temporairement le cordon avec le pointeur brut d'Enzyme
    // pour forcer la physique d'origine à lire le vecteur standard de la classe.
    this->phi_ptr = nullptr; 
    for(size_t i = 0; i < cell_nb; i++) {
        this->phi[i] = phi_val; 
    }
    
    // Appel de la méthode de simulation classique
    this->execute(); 
    
    // Calcul du coût sur l'état final
    double J_orig = 0.0;
    if (!obs_frames_storage.empty()) {
        J_orig = (obs_type == "gauges")
            ? this->compute_cost_gauges_ptr(this->obs_frames_storage.back().data())
            : this->compute_cost_frame_ptr(this->obs_frames_storage.back().data());
    }


    // --- 5. AFFICHAGE DU RÉSULTAT ---
    std::cout << "Phi = " << phi_val 
              << " | J_Enzyme = " << J_E 
              << " | J_Original = " << J_orig << std::endl;

    // Restauration de l'état initial du pointeur pour la suite du programme
    this->phi_ptr = old_phi_ptr;
}