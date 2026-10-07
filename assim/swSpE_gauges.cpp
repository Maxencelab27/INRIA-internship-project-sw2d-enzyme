/**
 * @file porosity_gauges.cpp
 * @brief Assimilation aux SONDES ponctuelles (obs_type=="gauges").
 *
 * Chemin d'observation SEPARE et ADDITIF : quand obs_type=="cells" (defaut), rien ici
 * n'est appele et l'assimilation SP/SP historique est strictement inchangee.
 *
 * Idee (approche 1) : la reference est une serie temporelle aux sondes physiques
 * (typiquement produite par le forward du modele fin SWES, meme probes_location.txt),
 * et l'operateur d'observation H_SP prend la valeur de la cellule SP contenant la sonde
 * (ordre 0, on reutilise probe_cell_idx deja rempli par swSpE::init). Aucune projection
 * de maillage, aucun facteur 1/phi : une sonde mesure la profondeur locale = h_SP.
 *
 * Le cout (compute_cost_gauges_ptr) est appele DANS la tape Enzyme (ckpt_segment_E) :
 * gauge_cell_idx est un index constant (enzyme_const), water_depth est actif -> Enzyme
 * differencie proprement par rapport a l'etat SP, puis a phi.
 *
 * Added for automatic differentiation; not part of native SW2D.
 */

#include <cstddef>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "swSpE.h"
#include "../../../common/param.h"
#include "../../../common/swLog.h"
#include "../../../common/geometry/swMesh.h"
#include "../../../common/geometry/swCell.h"


// ---------------------------------------------------------------------------
// locate_gauges : localise chaque sonde (probe_x/probe_y) dans le maillage SP.
// TOLERANT : une sonde hors maillage (les maillages SP et SWES different, certaines
// jauges peuvent tomber hors de l'emprise SP ou dans un trou) est SAUTEE avec un
// avertissement, PAS une erreur fatale. On memorise l'index d'ORIGINE de chaque sonde
// retenue (gauge_src_col) pour lire la bonne colonne dans la reference.
// ---------------------------------------------------------------------------
bool swSpE::locate_gauges(swMesh* mesh)
{
    gauge_cell_idx.clear();
    gauge_src_col.clear();

    const size_t np = probe_x.size();
    if (np == 0) {
        logError() << "obs_type=gauges : aucune sonde (probe_locations.txt absent ou vide).";
        return false;
    }

    size_t skipped = 0;
    for (size_t i = 0; i < np; ++i) {
        swCell* c = mesh->getCellFromCoordinates(probe_x[i], probe_y[i]);
        if (c == nullptr) {
            logWarn() << "obs_type=gauges : sonde " << i
                         << " (" << probe_x[i] << ", " << probe_y[i]
                         << ") hors du maillage SP -> ignoree.";
            skipped++;
            continue;
        }
        gauge_cell_idx.push_back(c->getIdx());
        gauge_src_col.push_back(i);
    }

    n_gauges = gauge_cell_idx.size();
    if (n_gauges == 0) {
        logError() << "obs_type=gauges : aucune sonde dans le maillage SP (" << np << " testee(s)).";
        return false;
    }
    if (gauge_w.size() != n_gauges)
        gauge_w.assign(n_gauges, 1.0);    // poids uniformes par defaut

    logInfo() << "obs_type=gauges : " << n_gauges << "/" << np
              << " sonde(s) localisee(s) dans le maillage SP (" << skipped << " hors maillage).";
    return true;
}


// ---------------------------------------------------------------------------
// load_gauge_observations : parse un fichier ProbesTimeSeries.txt (reference SWES).
//
// Format attendu (identique a la sortie storeResultProbes de SW2D) :
//   - N lignes d'en-tete "nom  x  y"
//   - 1 ligne colonnes : "time  <nom>_zb  <nom>_h  <nom>_q  <nom>_r  ..."
//   - lignes de donnees : time  zb1 h1 q1 r1  zb2 h2 q2 r2  ...
// Pour la sonde g (0-indexee), h = colonne (2 + 4*g).  On n'assimile que h.
//
// Ne garde que les temps dans [t_state_assimilation, t_max]. Remplit observations_map
// EXACTEMENT comme load_all_observations (meme contrat) mais avec des trames de taille
// n_gauges au lieu de cell_nb*3 -> toute la plomberie aval de main_enzyme.cpp
// (obs_times_vec / obs_frames_storage / current_obs_*) est reutilisee sans modification.
// ---------------------------------------------------------------------------
bool swSpE::load_gauge_observations(const std::string& file)
{
    if (n_gauges == 0 || gauge_src_col.size() != n_gauges) {
        logError() << "obs_type=gauges : sondes non localisees avant le chargement "
                      "(locate_gauges doit etre appele dans init).";
        return false;
    }

    std::ifstream in(file);
    if (!in.is_open()) {
        logError() << "obs_type=gauges : impossible d'ouvrir la reference sondes : " << file;
        return false;
    }

    const double t_begin = t_state_assimilation;
    const double t_end   = t_max;
    const double eps     = 1e-6;

    observations_map.clear();

    std::string line;
    bool header_passed = false;
    int loaded = 0, skipped = 0;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        // On saute tout jusqu'a la ligne de colonnes commencant par "time".
        if (!header_passed) {
            if (line.rfind("time", 0) == 0) header_passed = true;
            continue;
        }

        std::stringstream ss(line);
        std::vector<double> cols;
        double v;
        while (ss >> v) cols.push_back(v);

        const double t_obs = cols.empty() ? 0.0 : cols[0];
        if (cols.empty() || t_obs < t_begin - eps || t_obs > t_end + eps) { skipped++; continue; }

        std::vector<double> frame(n_gauges, 0.0);
        bool ok = true;
        for (size_t g = 0; g < n_gauges; ++g) {
            // colonne h de la sonde d'origine gauge_src_col[g] : time + 4/sonde, h a l'offset 2
            const size_t hc = 2 + 4 * gauge_src_col[g];
            if (hc >= cols.size()) { ok = false; break; }
            frame[g] = cols[hc];
        }
        if (!ok) continue;

        observations_map[t_obs] = std::move(frame);   // meme contrat que load_all_observations
        loaded++;
    }
    in.close();

    if (loaded == 0) {
        logError() << "obs_type=gauges : aucune observation dans [" << t_begin << ", "
                   << t_end << "] depuis " << file << " (" << skipped << " hors fenetre).";
        return false;
    }
    logInfo() << "obs_type=gauges : " << loaded << " observation(s) sonde(s) chargee(s) depuis "
              << file << " (" << skipped << " hors fenetre).";
    return true;
}


// ---------------------------------------------------------------------------
// compute_cost_gauges_ptr : cout L2 aux sondes.  obs_g pointe sur n_gauges valeurs (h).
// APPELE DANS LA TAPE ENZYME. gauge_cell_idx[g] est un index constant (enzyme_const) ;
// water_depth est actif -> derivation propre. Peu d'entrees actives (n_gauges cellules).
// ---------------------------------------------------------------------------
double swSpE::compute_cost_gauges_ptr(const double* obs_g)
{
    double J = 0.0;
    for (size_t g = 0; g < n_gauges; ++g) {
        const size_t c = gauge_cell_idx[g];
        const double dh = water_depth[c] - obs_g[g];
        J += gauge_w[g] * dh * dh;
    }
    // Meme normalisation optionnelle que compute_cost_frame_ptr (cost_normalize) :
    // current_n_obs ne derive pas de phi -> coherent avec le gradient (Taylor reste a 1).
    if (cost_normalize && current_n_obs > 0)
        J /= (double)current_n_obs;
    return J;
}
