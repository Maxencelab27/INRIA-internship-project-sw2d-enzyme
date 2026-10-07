// =====================================================================
//  main_enzyme.cpp — pilote de la couche d'assimilation / AD (Enzyme),
//  DÉCOUPLÉ du SW2D natif.
//
//  Différences avec le natif :
//   - instancie swSpE (dérivée AD) directement, pas via le factory natif ;
//   - lit les clés AD via assimConfig (le natif ne les connaît pas) et fait
//     lire au parseur natif une COPIE FILTRÉE d'input.solv (clés AD retirées)
//     dans un dossier temporaire ;
//   - pousse les valeurs AD dans swSpE après le importAttributesFromSolverParam
//     natif (final -> non surchargeable).
// =====================================================================
#include <chrono>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <fstream>

#include "../src/common/swFile.h"
#include "../src/common/swIO.h"
#include "../src/common/swLog.h"
#include "../src/common/geometry/swMesh.h"
#include "../src/common/swMessage.h"
#include "../src/common/version.h"
#include "../src/sw2dSolver/common/swSolverParam.h"

#include "swSpE.h"
#include "assimConfig.h"

namespace fs = std::filesystem;

int main(int argc, char **argv)
{
    std::cout << " lancement de main_e (decouple) " << std::endl;

    if (argc < 2 || argc > 5) {
        std::cerr << "sw2d_enzyme usage:\n"
                  << "  sw2d_enzyme <folder> [mode] [phi_init.csv] [alpha_reg]\n"
                  << "    mode : (none/forward) | taylor | gen-ref | gen-ref-e | minimize\n";
        return EXIT_FAILURE;
    }
    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        std::cout << "sw2d_enzyme : moteur d'assimilation de porosite (Enzyme, decouple).\n";
        return EXIT_SUCCESS;
    }
    if (!strcmp(argv[1], "-v") || !strcmp(argv[1], "--version")) {
        std::cout << "SW2D version " << sw2d_version << std::endl;
        return EXIT_SUCCESS;
    }

    const std::string user_folder = argv[1];
    if (!fs::exists(user_folder)) {
        std::cerr << "Folder " << user_folder << " does not exist\n";
        return EXIT_FAILURE;
    }

    swIO *m_io = new swIO(user_folder);
    m_io->createEmptyFolders(swFile::log_solver_file);
    swMessage::outfile = new std::ofstream(m_io->absolutePath(swFile::log_solver_file));

    std::string command_line;
    for (int i = 0; i < argc; ++i) command_line += (i ? " " : "") + std::string(argv[i]);
    logInfo() << "Command line:" << command_line;
    logInfo() << "Running version : " << sw2d_version;

    // =====================================================================
    //  CONFIG DÉCOUPLÉE : lire les clés AD + fabriquer une copie filtrée
    //  d'input.solv dans un dossier temporaire, lue par le parseur NATIF.
    // =====================================================================
    AssimConfig assim_cfg;
    const fs::path tmp = fs::path(user_folder) / ".assim_tmp";
    std::error_code ec;
    fs::remove_all(tmp, ec);
    fs::create_directories(tmp / "input", ec);
    // Miroir de tout input/ (mesh, maps, extraction_time, probes...) par symlink.
    for (const auto& e : fs::directory_iterator(fs::path(user_folder) / "input")) {
        fs::create_symlink(fs::absolute(e.path()),
                           tmp / "input" / e.path().filename(), ec);
    }
    // input.solv filtré (clés AD retirées) écrit par-dessus le symlink.
    fs::remove(tmp / "input" / "input.solv", ec);
    const std::string user_solv = (fs::path(user_folder) / "input" / "input.solv").string();
    const std::string tmp_solv  = (tmp / "input" / "input.solv").string();
    if (!assim_cfg.load_and_filter(user_solv, tmp_solv)) {
        logError() << "Echec lecture/filtrage de " << user_solv;
        return EXIT_FAILURE;
    }

    logInfo() << "BEGIN INIT";
    // Le parseur natif lit la copie filtrée depuis le dossier temporaire.
    swIO *io_cfg = new swIO(tmp.string());
    swSolverParam *m_solver_param = new swSolverParam(io_cfg);
    if (m_solver_param->readParameters() != EXIT_SUCCESS) {
        logError() << "Error in m_solver_param->readParameters()";
        return EXIT_FAILURE;
    }
    fs::remove_all(tmp, ec);   // temp plus nécessaire (tout est parsé)
    logInfo() << "END INIT";

    logInfo() << "BEGIN MESH BUILDING";
    swMesh* m_mesh = new swMesh(m_io);
    if (m_mesh->build() != EXIT_SUCCESS) {
        logError() << "Error while building mesh"; return EXIT_FAILURE;
    }
    if (m_solver_param->getHasQsection() && m_mesh->readQSections() != EXIT_SUCCESS) {
        logError() << "Error while reading qsections"; return EXIT_FAILURE;
    }
    logInfo() << "END MESH BUILDING";

    // --- Instanciation directe de swSpE (remplace swParent::create) ---
    logInfo() << "BEGIN MODEL BUILDING (swSpE)";
    swSpE* sp_model = new swSpE();
    sp_model->setIOPointer(m_io);
    logInfo() << "END MODEL BUILDING";

    if (!m_solver_param->getHotstart()) m_io->createEmptyOutputFolder(true);

    sp_model->importAttributesFromSolverParam(m_solver_param); // attrs NATIFS
    assim_cfg.apply_to(sp_model);                              // clés AD -> swSpE

    logInfo() << "BEGIN MODEL INITIALIZATION";
    if (sp_model->init(m_mesh) != EXIT_SUCCESS) {
        logError() << "Error while initialising the model"; return EXIT_FAILURE;
    }
    // Pondérations du coût (faisait partie de l'init AD-modifiée du natif ; ici
    // appelees explicitement cote AD, cell_nb etant connu apres init).
    sp_model->resolve_cost_weight(sp_model->h_weight_spec, sp_model->h_weight_map, "h_weight");
    sp_model->resolve_cost_weight(sp_model->q_weight_spec, sp_model->q_weight_map, "q_weight");
    sp_model->init_primal_memory();
    sp_model->save_initial_state();
    logInfo() << "END MODEL INITIALIZATION";

    const std::string mode = (argc >= 3) ? std::string(argv[2]) : std::string("");
    const std::size_t checkpoint_stride = assim_cfg.checkpoint_stride;

    // Override optionnel de la carte de porosite (ref / ebauche).
    {
        const bool is_genref = (mode == "gen-ref" || mode == "gen-ref-e");
        std::string ov, role;
        if (is_genref)               { ov = assim_cfg.porosity_ref_map;     role = "ref"; }
        else if (mode != "minimize") { ov = assim_cfg.porosity_ebauche_map; role = "ebauche"; }
        if (!ov.empty()) {
            const std::string p = m_io->absolutePath(ov);
            std::vector<double> tmpv((std::size_t)sp_model->cell_nb, 0.0);
            if (!sp_model->load_phi_csv(p, tmpv.data())) {
                logError() << "Echec chargement carte porosite '" << ov << "'."; return EXIT_FAILURE;
            }
            for (int i = 0; i < sp_model->cell_nb; i++) sp_model->phi_data()[i] = tmpv[(std::size_t)i];
            sp_model->computeEdgeFactors();
            logInfo() << "Carte porosite (" << role << ") imposee depuis " << p;
        }
    }

    // ---- MODE gen-ref-e : references auto-coherentes via le forward _E ----
    if (mode == "gen-ref-e") {
        const int n = sp_model->cell_nb;
        std::vector<double> phi_ref(n);
        for (int i = 0; i < n; i++) phi_ref[i] = sp_model->phi_data()[i];
        const bool ok = sp_model->generate_references_E(phi_ref.data(), checkpoint_stride);
        if (!ok) { logError() << "Echec gen-ref-e."; return EXIT_FAILURE; }
        logExit() << " Done (gen-ref-e)"; swMessage::finalizeLog(); return 0;
    }

    // ---- MODE forward (defaut) : forward _E seul ----
    if (mode == "" || mode == "forward") {
        const int n = sp_model->cell_nb;
        std::vector<double> phi_fwd(n);
        for (int i = 0; i < n; i++) phi_fwd[i] = sp_model->phi_data()[i];
        const bool ok = sp_model->run_forward_E(phi_fwd.data(), checkpoint_stride);
        if (!ok) { logError() << "Echec forward _E."; return EXIT_FAILURE; }
        logExit() << " Done (forward)"; swMessage::finalizeLog(); return 0;
    }

    // ---- MODE gen-ref : references via le moteur classique ----
    if (mode == "gen-ref") {
        if (!sp_model->begin_reference_generation()) {
            logError() << "Echec init gen-ref."; return EXIT_FAILURE;
        }
    }

    // Forward classique (gen-ref et taylor ; inutile pour minimize)
    if (mode != "minimize") {
        logInfo() << "BEGIN SW2D MODEL COMPUTATION";
        if (sp_model->execute() != EXIT_SUCCESS) {
            logError() << "Error while running the model"; return EXIT_FAILURE;
        }
        logInfo() << "END SW2D MODEL COMPUTATION";
    }

    if (mode == "gen-ref") {
        if (!sp_model->finalize_reference_generation()) {
            logError() << "Echec gen-ref."; return EXIT_FAILURE;
        }
        sp_model->reset_to_start();
        logExit() << " Done (gen-ref)"; swMessage::finalizeLog(); return 0;
    }

    sp_model->reset_to_start();

    // ================= Assimilation (taylor / minimize) =================
    logInfo() << "BEGIN ADJOINT COMPUTATION (ENZYME)";

    bool obs_loaded;
    if (sp_model->obs_type == "gauges") {
        obs_loaded = sp_model->load_gauge_observations(m_io->absolutePath("reference/gauges.txt"));
    } else {
        obs_loaded = sp_model->load_all_observations();
    }
    if (!obs_loaded) { logError() << "Chargement des observations echoue."; return EXIT_FAILURE; }
    if (sp_model->observations_map.empty()) {
        logError() << "Map d'observations vide -> gradient nul."; return EXIT_FAILURE;
    }

    // Pointeurs stables pour Enzyme
    const std::size_t n_obs = sp_model->observations_map.size();
    sp_model->obs_times_vec.clear();
    sp_model->obs_frames_storage.clear();
    sp_model->obs_frames_ptrs.clear();
    sp_model->obs_frames_storage.reserve(n_obs);
    for (auto const& kv : sp_model->observations_map) {
        sp_model->obs_times_vec.push_back(kv.first);
        sp_model->obs_frames_storage.push_back(kv.second);
    }
    for (std::size_t i = 0; i < n_obs; i++)
        sp_model->obs_frames_ptrs.push_back(sp_model->obs_frames_storage[i].data());
    sp_model->current_obs_times  = sp_model->obs_times_vec.data();
    sp_model->current_obs_frames = sp_model->obs_frames_ptrs.data();
    sp_model->current_n_obs      = n_obs;

    const int n = sp_model->cell_nb;
    std::vector<double> phi_map(n);
    for (int i = 0; i < n; i++) phi_map[i] = sp_model->phi_data()[i];

    // ---- MODE minimize ----
    if (mode == "minimize") {
        std::string phi_path;
        if (argc >= 4)                              phi_path = std::string(argv[3]);
        else if (!assim_cfg.porosity_ebauche_map.empty())
            phi_path = m_io->absolutePath(assim_cfg.porosity_ebauche_map);
        else                                        phi_path = m_io->absolutePath("phi_init_ones.csv");

        if (!sp_model->load_phi_csv(phi_path, phi_map.data())) {
            logError() << "Echec chargement carte porosite initiale."; return EXIT_FAILURE;
        }
        double alpha_reg = assim_cfg.alpha_reg;
        if (argc >= 5) alpha_reg = std::atof(argv[4]);
        sp_model->run_minimization(phi_map.data(), assim_cfg.max_it, checkpoint_stride, alpha_reg);
        logExit() << " Done (minimize)"; swMessage::finalizeLog(); return 0;
    }

    // ---- MODE taylor : un gradient adjoint (+ Taylor), export CSV ----
    std::vector<double> dJ_dPhi(n, 0.0);
    sp_model->reset_to_start();
    logInfo() << "Calcul du gradient (checkpointe, stride=" << checkpoint_stride << ")...";
    sp_model->get_spatial_gradient_checkpointed(phi_map.data(), dJ_dPhi.data(), checkpoint_stride);

    std::ofstream out_csv(m_io->absolutePath("gradient_porosity.csv"));
    if (out_csv.is_open()) {
        out_csv << "cell_id,phi,gradient\n";
        for (int i = 0; i < n; i++) out_csv << i << "," << phi_map[i] << "," << dJ_dPhi[i] << "\n";
        out_csv.close();
        logInfo() << "Gradient exporte : gradient_porosity.csv";
    }

    logExit() << " Done";
    swMessage::finalizeLog();
    return 0;
}
