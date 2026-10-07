// ===========================================================================
//  porosity_minimize.cpp
//
//  Data-assimilation drivers for the single-porosity (SP) model:
//    1. export_reference_observations() — turns the forward output maps into
//       reference observations (reference/ResultHydro_p*.txt) sampled every
//       dt_map_ad inside the assimilation window [t_state_assimilation, t_max].
//    2. load_phi_csv()                  — reads a starting porosity map (CSV).
//    3. run_minimization()              — minimizes the L2 misfit cost J(phi)
//       with the PSQN dense quasi-Newton (BFGS) optimizer.
//
//  Cost / gradient consistency
//  ---------------------------
//  The optimizer needs f(phi) and ∇f(phi) to be consistent. We reuse the
//  checkpointed adjoint machinery that the Taylor test already validates:
//    - grad(phi): build_checkpoint_records(phi) (inside
//                 get_spatial_gradient_checkpointed) records a FIXED
//                 trajectory (dt/BC/obs) anchored at phi, then Enzyme returns
//                 ∂J/∂phi for that frozen trajectory.
//    - func(phi): replay_trajectory_cost(phi) evaluates J on the SAME frozen
//                 trajectory → smooth in phi (no CFL kinks) and exactly the
//                 function whose gradient grad() returns.
//  PSQN::bfgs always calls grad() at an iterate before func() probes the line,
//  so the records are anchored before any func() evaluation.
// ===========================================================================

#include "swSpE.h"
#include "../swParent.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstdio>
#include <limits>
#include <algorithm>
#include <cmath>
#include <map>

#include "../../../common/swLog.h"
#include "../../../common/swIO.h"
#include "../../../common/swFile.h"
#include "../../../common/swString.h"

#include "psqn-bfgs.h"

// L-BFGS-B a contraintes de boite (LBFGSpp, header-only, backend Eigen). Utilise quand
// minimize_algo="lbfgsb" : impose nativement phi in [phi_lower_bound, phi_upper_bound]
// (pas de transform sigmoide, donc pas de saturation dphi/dm). Cf recap LaTeX livre a part.
#include <Eigen/Core>
#include <LBFGSB.h>

// NB : l'optimiseur Adam est implemente A LA MAIN (boucle projetee) dans
// run_minimization, et NON via optim::gd : optim impose les bornes par un transform
// interne (sigmoide) qui reintroduirait la saturation dphi/dm->0. La lib optim reste
// liee (CMake) mais n'est plus utilisee ici.

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// 1. Reference observation generation
//
//    Les references ne sont PAS des copies des cartes output : elles sont
//    ecrites directement par le solveur pendant la passe forward gen-ref, a
//    chaque dt_map_ad de la fenetre [t_state_assimilation, t_max] (cf.
//    storeReferenceMap, appelee depuis execute()). Les cartes output restent,
//    elles, sur la grille dt_map. dt_map_ad < dt_map est autorise.
// ---------------------------------------------------------------------------
bool swSpE::begin_reference_generation()
{
    const double dt_ref = dt_map_ad;            // already resolved (== dt_map if unset)

    if (dt_ref <= 0.) {
        logError() << "begin_reference_generation : dt_map_ad invalide (" << dt_ref << ")";
        return false;
    }

    // Le solveur cale ses pas de temps (et donc les references) sur les multiples
    // de dt_map_ad depuis t=0. t_state_assimilation doit etre sur cette grille,
    // sinon aucune reference ne tombe exactement sur la borne de la fenetre.
    double ratio = t_state_assimilation / dt_ref;
    if (std::abs(ratio - std::round(ratio)) > 1e-6) {
        logWarn() << "begin_reference_generation : t_state_assimilation (" << t_state_assimilation
                  << ") n'est pas un multiple de dt_map_ad (" << dt_ref << ") ;"
                  << " la premiere reference sera prise au premier multiple de dt_map_ad dans la fenetre.";
    }

    // (Re)cree le dossier reference/
    const std::string ref_dir = m_io->absolutePath("reference");
    std::error_code ec;
    fs::remove_all(ref_dir, ec);
    fs::create_directories(ref_dir, ec);
    if (ec) {
        logError() << "begin_reference_generation : impossible de creer " << ref_dir
                   << " (" << ec.message() << ")";
        return false;
    }

    n_ref_written = 0;
    gen_ref_mode  = true;

    logInfo() << "Generation des references dans la fenetre [" << t_state_assimilation
              << "s, " << t_max << "s] tous les " << dt_ref << "s (dt_map=" << dt_map
              << "s pour les cartes output) -> " << ref_dir;
    return true;
}

bool swSpE::finalize_reference_generation()
{
    gen_ref_mode = false;
    logInfo() << n_ref_written << " observation(s) de reference generee(s) directement"
              << " dans reference/ (tous les " << dt_map_ad << "s).";
    if (n_ref_written == 0) {
        logError() << "Aucune reference generee : verifier la fenetre [" << t_state_assimilation
                   << "s, " << t_max << "s] et dt_map_ad (" << dt_map_ad << ").";
        return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// 2. Starting porosity map : CSV ("cell_id,phi") OU format SP_cell_porosity_map.txt
//    (Unif/Distrib). Auto-detection : un fichier commencant par "Unif" est lu comme SP
//    (delegue a swIO::readVectorParameterFile), sinon comme CSV. Permet de re-injecter
//    directement une carte optimisee (ecrite au format SP) comme ebauche/ref.
// ---------------------------------------------------------------------------
bool swSpE::load_phi_csv(const std::string& filepath, double* phi_out)
{
    std::ifstream file(filepath);
    if (!file.is_open()) {
        logError() << "load_phi_csv : impossible d'ouvrir " << filepath;
        return false;
    }

    // --- Detection du format : SP Unif/Distrib (1ere ligne non vide = "Unif ...") ---
    {
        std::string first;
        while (std::getline(file, first)) {
            size_t p = first.find_first_not_of(" \t\r\n");
            if (p != std::string::npos) { first = first.substr(p); break; }
        }
        if (first.rfind("Unif", 0) == 0 || first.rfind("unif", 0) == 0) {
            file.close();
            // Format SP_cell_porosity_map.txt : lecteur canonique (gere Unif 0 et 1).
            std::vector<std::string> names;
            std::vector<std::vector<double>> arr;
            swIO::readVectorParameterFile(filepath, 1, cell_nb, names, arr);
            if (arr.empty() || arr[0].size() != cell_nb) {
                logError() << "load_phi_csv (format SP) : nb de valeurs != cell_nb dans " << filepath;
                return false;
            }
            for (size_t i = 0; i < cell_nb; i++) phi_out[i] = arr[0][i];
            logInfo() << "Carte de porosite chargee (format SP Unif/Distrib) depuis " << filepath
                      << " (" << cell_nb << " cellules).";
            return true;
        }
        file.clear();
        file.seekg(0);   // rembobiner pour la lecture CSV ci-dessous
    }

    std::string line;
    // En-tete optionnel : si la premiere ligne n'est pas numerique, on la saute.
    std::streampos data_start = file.tellg();
    if (std::getline(file, line)) {
        std::string trimmed = line;
        bool looks_numeric = !trimmed.empty() &&
            (std::isdigit((unsigned char)trimmed[0]) || trimmed[0] == '-' || trimmed[0] == '+');
        if (looks_numeric) file.seekg(data_start);   // pas d'en-tete : rembobiner
    }

    size_t count = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::stringstream ss(line);
        long id;
        double val;
        if (!(ss >> id >> val)) continue;
        if (id < 0 || (size_t)id >= cell_nb) {
            logError() << "load_phi_csv : cell_id hors limites (" << id << ") dans " << filepath;
            return false;
        }
        phi_out[id] = val;
        count++;
    }
    file.close();

    if (count != cell_nb) {
        logError() << "load_phi_csv : " << count << " valeur(s) lue(s) mais "
                   << cell_nb << " cellule(s) attendue(s) (" << filepath << ")";
        return false;
    }
    logInfo() << "Carte de porosite initiale chargee depuis " << filepath
              << " (" << count << " cellules).";
    return true;
}

// ---------------------------------------------------------------------------
// 3. Minimization driver (PSQN dense BFGS)
// ---------------------------------------------------------------------------

// Console reporter: prints one line per BFGS iteration.
namespace {
struct MinimizeReporter : public PSQN::dummy_reporter {
    static void line_search
    (int const /*trace*/, PSQN::psqn_uint const iteration,
     PSQN::psqn_uint const n_eval, PSQN::psqn_uint const n_grad,
     double const fval_old, double const fval, bool const successful,
     double const step_size, double const * /*new_x*/,
     PSQN::psqn_uint const /*n_global*/)
    {
        std::cout << "[BFGS] it=" << std::setw(3) << iteration
                  << "  J=" << std::setprecision(10) << std::scientific << fval
                  << "  (prev " << fval_old << ")"
                  << "  step=" << step_size
                  << "  n_eval=" << n_eval << "  n_grad=" << n_grad
                  << (successful ? "" : "   [line-search FAILED]")
                  << std::endl;
    }
};
} // namespace

// Local problem wrapper handed to PSQN::bfgs.
//
// Bound handling (changement de variable m -> phi, pilote par phi_transform):
// PSQN::bfgs est NON CONTRAINT ; la variable d'optimisation est m (la porosite
// brute, non bornee). Le passage m -> phi se fait par swSpE::phi_correct :
//     none      : phi = m
//     clamp     : phi = clamp(m, 0, 1.05)
//     tanh      : phi = 1/(1+exp(-m))
//     normalize : phi = (m-lo)/(hi-lo), lo/hi = min/max des cellules libres (par iteree)
// appliquE DANS la fonction differenciee par Enzyme (ckpt_segment_E /
// build_checkpoint_records). Donc Enzyme rend directement dJ_data/dm (le terme
// dphi/dm est dans la tape). Seul le ridge ajoute sa derivee via dphi_dm() (chain
// rule analytique, hors tape).
namespace {

// Ecrit une carte de porosite EXACTEMENT au format SP_cell_porosity_map.txt
// (Unif/Distrib, 1 colonne) -> rechargeable tel quel par swIO::readVectorParameterFile
// (cp vers input/SP_cell_porosity_map.txt). Structure byte-identique a l'entree :
//
//   Unif<TAB>0          <- mode distribue (n valeurs)
//   ==== Default Param
//   Phi
//   1
//   ==== Distrib
//   <phi_0> ... <phi_{n-1}>   <- n valeurs (1 par ligne)
//
// NB : le lecteur lit des POSITIONS FIXES (readLine ne saute aucun commentaire) -> ne
// rien ajouter/retirer. Aucun contexte/commentaire inseré (format strict).
static void write_sp_porosity_map(const std::string& path, const double* phi, size_t n)
{
    std::ofstream o(path);
    if (!o.is_open()) return;
    o << "Unif\t0\n";
    o << "==== Default Param\n";
    o << "Phi\n";
    o << "1\n";
    o << "==== Distrib\n";
    o << std::setprecision(10);
    for (size_t i = 0; i < n; i++) o << phi[i] << "\n";
}

struct PorosityProblem : public PSQN::problem {
    swSpE*  solver;
    size_t n;                      // total number of cells
    size_t stride;
    double alpha;                  // Tikhonov ridge weight
    const double* phi_b;           // background (= initial map), size n
    const std::vector<size_t>& free_idx;  // indices of the cells actually optimized
    std::vector<double> m_buf;     // full OPTIMIZATION variable m (size n). phi = solver->phi_correct(m)
    std::vector<double> gm_buf;    // full dJ_data/dm from the adjoint (size n)

    // --- ZONATION (optionnelle) : groupe les cellules libres par valeur d'ebauche (phi_b)
    // identique. La variable d'optimisation devient 1 valeur PAR ZONE (au lieu d'1 par cellule
    // libre) : une region uniforme de l'ebauche bouge en bloc. Le gradient d'une zone = SOMME
    // des gradients de ses cellules -> l'interieur masque herite du gradient fort du bord.
    bool zonation_on = false;
    std::vector<size_t> zone_of_free;   // [j] -> index de zone de la j-eme cellule libre
    std::vector<double> zone_value;     // [z] -> valeur d'ebauche representative de la zone z
    size_t n_zones = 0;
    std::vector<double> gcell_tmp;      // gradient par cellule libre (avant reduction par zone)

    // Suivi de convergence : J au point de depart (1ere eval de grad) et compteur d'eval.
    double J_init   = 0.0;
    bool   J_init_set = false;
    int    grad_calls = 0;

    PorosityProblem(swSpE* s, size_t st, double a, const double* pb,
                    const std::vector<size_t>& fidx, bool zon = false)
        : solver(s), n(s->cell_nb), stride(st),
          alpha(a), phi_b(pb), free_idx(fidx),
          m_buf(s->cell_nb), gm_buf(s->cell_nb, 0.0),
          zonation_on(zon), gcell_tmp(fidx.size(), 0.0) {
        if (zonation_on) build_zones();
    }

    // Groupe les cellules libres par valeur d'ebauche (phi_b) identique -> zones.
    void build_zones() {
        zone_of_free.assign(free_idx.size(), 0);
        zone_value.clear();
        std::map<long, size_t> val2zone;               // cle = valeur quantifiee -> index de zone
        for (size_t j = 0; j < free_idx.size(); j++) {
            double v = phi_b[free_idx[j]];
            long key = (long)std::llround(v * 1e6);     // quantifie pour absorber le bruit flottant
            auto it = val2zone.find(key);
            size_t z;
            if (it == val2zone.end()) { z = zone_value.size(); val2zone[key] = z; zone_value.push_back(v); }
            else                        z = it->second;
            zone_of_free[j] = z;
        }
        n_zones = zone_value.size();
    }

    // The optimization runs over the FREE cells only (mask), ou par ZONE si zonation.
    PSQN::psqn_uint size() const override {
        return (PSQN::psqn_uint)(zonation_on ? n_zones : free_idx.size());
    }

    // theta = la variable d'optimisation m sur les cellules libres.
    // -> m_buf complet : cellules fixees = phi_b, cellules libres = theta. Le passage
    //    m -> phi est fait DANS la tape Enzyme (swSpE::phi_correct appele dans
    //    ckpt_segment_E / build_checkpoint_records), donc on passe m_buf "tel quel"
    //    comme entree differenciee.
    void map_theta_to_phi(double const* theta) {
        for (size_t i = 0; i < n; i++) m_buf[i] = phi_b[i];          // fixed cells = background
        if (zonation_on) {
            for (size_t j = 0; j < free_idx.size(); j++)
                m_buf[free_idx[j]] = theta[zone_of_free[j]];        // free cells = valeur de zone
        } else {
            for (size_t j = 0; j < free_idx.size(); j++)
                m_buf[free_idx[j]] = theta[j];                      // free cells = m
        }
    }

    // Ridge penalty (alpha/2) * sum (phi - phi_b)^2, avec phi = phi_correct(m).
    double ridge_penalty() const {
        if (alpha <= 0.0) return 0.0;
        double s = 0.0;
        for (size_t j = 0; j < free_idx.size(); j++) {
            double phi = solver->phi_correct(m_buf[free_idx[j]]);
            double d = phi - phi_b[free_idx[j]]; s += d * d;
        }
        return 0.5 * alpha * s;
    }

    // Penalite L1 de parcimonie : alpha_l1 * sum |phi - phi_b|. Repousse les nombreuses
    // petites deviations vers phi_b (nettoie la baisse diffuse) en autorisant quelques
    // grosses (le barrage). Derivee = alpha_l1 * sign(phi-phi_b) * dphi/dm (ajoutee dans grad).
    double l1_penalty() const {
        const double a_l1 = solver->get_alpha_l1();
        if (a_l1 <= 0.0) return 0.0;
        double s = 0.0;
        for (size_t j = 0; j < free_idx.size(); j++) {
            double phi = solver->phi_correct(m_buf[free_idx[j]]);
            s += std::fabs(phi - phi_b[free_idx[j]]);
        }
        return a_l1 * s;
    }

    // Regularisation SPATIALE lisse (H1 / Laplacien de graphe) :
    //   0.5 * alpha_spatial * sum_{aretes internes} (phi_L - phi_R)^2.
    // Somme sur TOUTES les aretes internes (couple cellules libres ET fixes) -> propage la valeur
    // de bord OBSERVEE vers l'interieur hydrauliquement masque du bloc. Lisse (differentiable
    // partout) -> compatible line-search L-BFGS-B (contrairement a TV/L1 non-lisses).
    double spatial_penalty() const {
        const double a = solver->get_alpha_spatial();
        if (a <= 0.0) return 0.0;
        const std::vector<size_t>& cL = solver->edge_int_celL_idx;
        const std::vector<size_t>& cR = solver->edge_int_celR_idx;
        double s = 0.0;
        for (size_t e = 0; e < cL.size(); e++) {
            double phiL = solver->phi_correct(m_buf[cL[e]]);
            double phiR = solver->phi_correct(m_buf[cR[e]]);
            double d = phiL - phiR; s += d * d;
        }
        return 0.5 * a * s;
    }

    // CHEAT/diagnostic : terme d'identifiabilite (alpha/2) * sum (phi - phi_ref)^2 ancre sur
    // la carte VERITE (porosity_ref_map). alpha = input.solv 'alpha_phi_ref'. Rend le probleme
    // bien pose (le minimum du cout est la verite) -> sert a tester la machinerie d'optim.
    double phi_ref_penalty() const {
        const double a = solver->get_alpha_phi_ref();
        if (a <= 0.0 || phi_ref.size() != n) return 0.0;
        double s = 0.0;
        for (size_t j = 0; j < free_idx.size(); j++) {
            double phi = solver->phi_correct(m_buf[free_idx[j]]);
            double d = phi - phi_ref[free_idx[j]]; s += d * d;
        }
        return 0.5 * a * s;
    }

    // f(theta) = J_data + ridge + L1 + spatial(H1) + (cheat) phi_ref, sur la trajectoire figee (records ancres par le dernier grad()).
    double func(double const* theta) override {
        map_theta_to_phi(theta);
        return solver->replay_trajectory_cost(m_buf.data()) + ridge_penalty() + l1_penalty() + spatial_penalty() + phi_ref_penalty();
    }

    // f(theta) et ∇_m f sur les cellules libres. ENZYME rend deja dJ_data/dm
    // (le changement de variable phi_correct est DANS la tape) -> AUCUN produit par
    // dphi/dtheta ici pour le terme donnees. Le ridge ajoute sa derivee analytique
    // (chain via dphi_dm() du mode actif).
    double grad(double const* theta, double* gr) override {
        map_theta_to_phi(theta);
        solver->get_spatial_gradient_checkpointed(m_buf.data(), gm_buf.data(), stride);
        double J_data = solver->replay_trajectory_cost(m_buf.data());

        // Boucle unique : gradient (gr) + metriques de suivi (phi et gradient-m sur les libres).
        const double a_l1 = solver->get_alpha_l1();
        // CHEAT/diagnostic : gradient du terme (alpha/2)*||phi - phi_ref||^2 = alpha*(phi-phi_ref)*dphi/dm.
        const double a_ref = solver->get_alpha_phi_ref();
        const bool   ref_ok = (a_ref > 0.0 && phi_ref.size() == n);
        // Gradient de la regularisation spatiale (H1) : d(0.5*a*sum(phi_L-phi_R)^2)/dphi_i
        //   = a * sum_{aretes touchant i} (phi_i - phi_voisin). Pre-calcule par cellule (sur
        //   TOUTES les aretes, y compris vers des cellules fixes) puis chaine par dphi/dm ci-dessous.
        const double a_sp = solver->get_alpha_spatial();
        std::vector<double> gspat;
        if (a_sp > 0.0) {
            gspat.assign(n, 0.0);
            const std::vector<size_t>& cL = solver->edge_int_celL_idx;
            const std::vector<size_t>& cR = solver->edge_int_celR_idx;
            for (size_t e = 0; e < cL.size(); e++) {
                double phiL = solver->phi_correct(m_buf[cL[e]]);
                double phiR = solver->phi_correct(m_buf[cR[e]]);
                double d = a_sp * (phiL - phiR);
                gspat[cL[e]] += d; gspat[cR[e]] -= d;
            }
        }
        double pmin = 2.0, pmax = -1.0, phi_sum = 0.0;
        double gm_l2 = 0.0, gm_linf = 0.0; long gm_argmax = -1;
        long n_lt_0p5 = 0, n_lt_0p05 = 0, n_upper = 0;
        for (size_t j = 0; j < free_idx.size(); j++) {
            size_t i   = free_idx[j];
            double m   = m_buf[i];
            double phi = solver->phi_correct(m);
            // dphi/dm analytique du mode actif (none/clamp/tanh/normalize). Le terme
            // dJ_data/dm est deja rendu par Enzyme (phi_correct est DANS la tape) ; ici
            // on ajoute les derivees analytiques du ridge (L2) et de la parcimonie (L1).
            double dphidm = solver->dphi_dm(m);
            double d_l1 = 0.0;                      // derivee L1 : alpha_l1 * sign(phi-phi_b) * dphi/dm
            if (a_l1 > 0.0) {
                double diff = phi - phi_b[i];
                if (diff >  1e-12) d_l1 =  a_l1 * dphidm;
                else if (diff < -1e-12) d_l1 = -a_l1 * dphidm;
            }
            // Gradient PAR CELLULE LIBRE (tampon) : data + dridge + dL1 + dspatial + d(cheat).
            double gval = gm_buf[i] + alpha * (phi - phi_b[i]) * dphidm + d_l1;
            if (a_sp > 0.0) gval += gspat[i] * dphidm;                    // + d(spatial H1)/dm
            if (ref_ok) gval += a_ref * (phi - phi_ref[i]) * dphidm;      // + d(cheat phi_ref)/dm
            gcell_tmp[j] = gval;
            if (phi < pmin) pmin = phi;
            if (phi > pmax) pmax = phi;
            phi_sum += phi;
            if (phi < 0.5)   n_lt_0p5++;
            if (phi < 0.05)  n_lt_0p05++;
            if (phi > 1.049) n_upper++;                // colle a la borne sup du clamp (1.05)
            double ag = std::fabs(gval);
            gm_l2 += gval * gval;
            if (ag > gm_linf) { gm_linf = ag; gm_argmax = (long)i; }
        }
        gm_l2 = std::sqrt(gm_l2);

        // Reduction : gradient par cellule libre -> variable d'optimisation (gr).
        // Sans zonation : gr[j] = gcell_tmp[j]. Avec zonation : gr[z] = SOMME des cellules
        // de la zone z (le bloc entier suit le gradient agrege, domine par le bord observable).
        if (zonation_on) {
            for (size_t z = 0; z < n_zones; z++) gr[z] = 0.0;
            for (size_t j = 0; j < free_idx.size(); j++) gr[zone_of_free[j]] += gcell_tmp[j];
        } else {
            for (size_t j = 0; j < free_idx.size(); j++) gr[j] = gcell_tmp[j];
        }
        const double phi_mean = free_idx.empty() ? 0.0 : phi_sum / (double)free_idx.size();
        const double ridge = ridge_penalty();
        double J = J_data + ridge + l1_penalty() + phi_ref_penalty();

        // --- Suivi de convergence : J_init (1ere eval) + differences relatives signees ---
        if (!J_init_set) { J_init = J; J_init_set = true; }
        ++grad_calls;
        const double rel      = (J_init != 0.0) ? (J - J_init) / J_init * 100.0 : 0.0;
        const double rel_prev = (J_prev_set && J_prev != 0.0) ? (J - J_prev) / J_prev * 100.0 : 0.0;
        const char rel_sign = (rel >= 0.0) ? '+' : '-';
        // Ecart parametrique a la carte VERITE (porosity_ref_map), affiche en direct
        // a chaque iteree (== colonne sse_phi_ref du descent_log.csv). -1 => non dispo.
        const double sse_ref = sse_to_ref();
        std::cout << "[Minim] grad #" << grad_calls
                  << " : J_current = " << std::scientific << std::setprecision(6) << J
                  << " | J_init = " << J_init
                  << " | rel = " << rel_sign << std::fixed << std::setprecision(2) << std::fabs(rel) << " %"
                  << std::scientific
                  << " | phi_free in [" << pmin << ", " << pmax << "]";
        if (sse_ref >= 0.0)
            std::cout << " | sse_phi_ref = " << std::setprecision(6) << sse_ref;
        std::cout << std::endl;

        // --- Fichier de suivi de descente (CSV) + carte phi de l'iteree ---
        write_descent_row(grad_calls, J, J_data, ridge, rel, rel_prev,
                          gm_l2, gm_linf, gm_argmax,
                          pmin, pmax, phi_mean, n_lt_0p5, n_lt_0p05, n_upper);
        dump_phi_map(grad_calls);

        J_prev = J; J_prev_set = true;

        // Apres la 1ere evaluation, on coupe les bannieres verbeuses (computeEdgeFactors,
        // [Checkpoint], [Window], GRADIENT CHECKPOINTED) : chaque iteree suivante reste terse.
        solver->set_quiet_logs(true);
        return J;
    }

    // ------------------------------------------------------------------
    // Suivi de descente : un CSV (1 ligne/eval) + 1 carte phi/eval.
    // descent_dir est cree par run_minimization. Ouverture paresseuse du CSV.
    // ------------------------------------------------------------------
    std::string   descent_dir;
    std::ofstream descent_csv;
    double J_prev = 0.0; bool J_prev_set = false;
    // Carte VERITE (porosity_ref_map) ; vide => colonne sse_phi_ref non remplie.
    // sse_phi_ref = sum_i (phi_i - phi_ref_i)^2 : ecart parametrique a la verite.
    std::vector<double> phi_ref;

    // phi physique courant complet : cellules fixees = phi_b, libres = phi_correct(m).
    void current_phi(std::vector<double>& phi_out) const {
        phi_out.assign(phi_b, phi_b + n);
        for (size_t j = 0; j < free_idx.size(); j++)
            phi_out[free_idx[j]] = solver->phi_correct(m_buf[free_idx[j]]);
    }
    // Somme des carres des ecarts a la carte verite (sur TOUTES les cellules).
    double sse_to_ref() const {
        if (phi_ref.size() != n) return -1.0;   // -1 => non disponible
        std::vector<double> phi; current_phi(phi);
        double s = 0.0;
        for (size_t i = 0; i < n; i++) { double d = phi[i] - phi_ref[i]; s += d * d; }
        return s;
    }

    void write_descent_row(int it, double J, double Jdata, double ridge,
                           double rel_init, double rel_prev,
                           double gm_l2, double gm_linf, long gm_argmax,
                           double phi_min, double phi_max, double phi_mean,
                           long n_lt_0p5, long n_lt_0p05, long n_upper)
    {
        if (descent_dir.empty()) return;
        if (!descent_csv.is_open()) {
            descent_csv.open(descent_dir + "/descent_log.csv");
            descent_csv << "iter,J_total,J_data,ridge,rel_init_pct,rel_prev_pct,"
                           "gradJm_L2,gradJm_Linf,gradJm_argmax,"
                           "gradJphys_L2,gradJphys_max,gradJphys_argmax,dh_peak,"
                           "phi_min,phi_max,phi_mean,n_phi_lt_0p5,n_phi_lt_0p05,"
                           "n_phi_upper_clamp,adam_step,sse_phi_ref\n";
        }
        descent_csv << it << ","
                    << std::scientific << std::setprecision(10)
                    << J << "," << Jdata << "," << ridge << ","
                    << std::setprecision(6) << rel_init << "," << rel_prev << ","
                    << gm_l2 << "," << gm_linf << "," << gm_argmax << ","
                    << solver->get_diag_grad_l2() << "," << solver->get_diag_grad_max() << ","
                    << solver->get_diag_grad_argmax() << "," << solver->get_diag_dh_peak() << ","
                    << phi_min << "," << phi_max << "," << phi_mean << ","
                    << n_lt_0p5 << "," << n_lt_0p05 << "," << n_upper << ","
                    << solver->get_adam_step_size() << ","
                    << std::scientific << std::setprecision(10) << sse_to_ref() << "\n";
        descent_csv.flush();   // pour suivre en direct pendant le run
    }

    void dump_phi_map(int it)
    {
        if (descent_dir.empty()) return;
        char name[64];
        // Format SP_cell_porosity_map.txt (Unif/Distrib) -> directement re-utilisable
        // comme carte de porosite (et traçable par les outils mesh_plot).
        std::snprintf(name, sizeof(name), "/phi_iter/phi_iter_%04d.txt", it);
        // phi physique : cellules fixees = phi_b, cellules libres = phi_correct(m).
        std::vector<double> phi_out(phi_b, phi_b + n);
        for (size_t j = 0; j < free_idx.size(); j++)
            phi_out[free_idx[j]] = solver->phi_correct(m_buf[free_idx[j]]);
        write_sp_porosity_map(descent_dir + name, phi_out.data(), n);
    }
};
} // namespace

double swSpE::run_minimization(double* phi_start, size_t max_it, size_t checkpoint_stride,
                              double alpha_reg)
{
    // Le test de Taylor (~1 forward par segment) serait relance a chaque gradient :
    // on le desactive pendant la minimisation.
    taylor_test_after_gradient = false;

    // Affichage : 1ere evaluation de gradient verbeuse (setup complet), puis terse.
    // Remis a false ici pour que le forward ebauche + la 1ere iteree s'affichent.
    minim_quiet_logs = false;

    // Changement de variable m -> phi applique DANS la tape Enzyme via swSpE::phi_correct,
    // pilote par 'phi_transform' (none/clamp/tanh/normalize). L'optimiseur travaille sur m
    // (non borne) et Enzyme rend directement dJ_data/dm. Pour 'normalize', les bornes lo/hi
    // (min/max des cellules libres) sont recalculees et figees a chaque iteree dans
    // get_spatial_gradient_checkpointed (update_phi_norm_bounds), avant la tape.

    // Background phi_b = carte de depart (figee). Le terme ridge tire les cellules
    // non informees vers phi_b et rend le probleme bien pose (Hessien + alpha*I).
    std::vector<double> phi_background(phi_start, phi_start + cell_nb);

    // --- Masque d'identifiabilite ---
    // Forward ebauche avec la carte de depart (ici phi=1 partout). Les cellules
    // encore SECHES a la fin (h <= h_dry) ne sont jamais atteintes par l'eau, meme
    // en flux libre (elles sont seches a cause de leur topographie z, hors riviere) :
    // leur porosite n'a aucun effet observable -> on les FIXE a phi_b. Les cellules
    // mouillees (dont l'emplacement du barrage, mouille sous phi=1) restent libres.
    const double h_dry = 1e-4;
    for (size_t i = 0; i < cell_nb; i++) phi[i] = phi_start[i];
    phi_ptr = phi.data();
    // FIX A : recalcul des facteurs d'arête (1/phi) apres avoir change phi vers la
    // carte de depart, sinon le forward ebauche utilise les facteurs figes a l'init
    // (carte d'entree) -> incoherence -> blow-up (ex: crash a t~2800s observe).
    computeEdgeFactors();
    reset_to_start();

    std::vector<size_t> free_idx;
    free_idx.reserve(cell_nb);
    if (mask_dry_cells) {
        // Masque d'identifiabilite ACTIF : on fige les cellules seches (h<=h_dry) a la fin
        // du forward ebauche (leur porosite n'a aucun effet observable).
        logInfo() << "Masque d'identifiabilite ACTIF (mask_dry_cells=1) : forward ebauche...";
        if (execute() != EXIT_SUCCESS) {
            logError() << "run_minimization : echec du forward ebauche (masque).";
        }
        size_t n_dry = 0;
        for (size_t i = 0; i < cell_nb; i++) {
            if (water_depth[i] > h_dry) free_idx.push_back(i);
            else n_dry++;
        }
        reset_to_start();
        logInfo() << "Masque : " << n_dry << " cellules seches FIXEES a phi_b, "
                  << free_idx.size() << " cellules libres optimisees (sur " << cell_nb << ").";
    } else {
        // Defaut : pas de masque -> TOUTES les cellules sont optimisees (seches incluses).
        // Les cellules a gradient nul (seches) ne bougeront simplement pas (Adam : pas nul).
        for (size_t i = 0; i < cell_nb; i++) free_idx.push_back(i);
        logInfo() << "Masque d'identifiabilite DESACTIVE (mask_dry_cells=0) : toutes les "
                  << cell_nb << " cellules optimisees (cellules seches incluses).";
    }

    // --- ZONATION : restriction spatiale des cellules libres (porosity_free_zone_map) ---
    // Carte CSV (cell_id,value) : cellule LIBRE ssi value < 0.5 ; hors zone => FIXEE a phi_b.
    // Combine avec le masque ci-dessus (intersection). Combat l'etalement diffus de la baisse
    // de porosite : l'optimiseur ne peut agir QUE dans la zone (ex. l'emprise du barrage).
    if (!porosity_free_zone_map_path.empty()) {
        std::vector<double> zone(cell_nb, 1.0);
        const std::string zpath = m_io->absolutePath(porosity_free_zone_map_path);
        if (!load_phi_csv(zpath, zone.data())) {
            logError() << "run_minimization : zone map illisible (" << zpath
                       << ") -> zonation ignoree, toutes cellules du masque conservees.";
        } else {
            std::vector<size_t> zoned;
            zoned.reserve(free_idx.size());
            for (size_t idx : free_idx)
                if (zone[idx] < 0.5) zoned.push_back(idx);   // free iff dans la zone
            logInfo() << "ZONATION active (porosity_free_zone_map=" << zpath << ") : "
                      << zoned.size() << " cellules libres dans la zone / "
                      << free_idx.size() << " (avant zonation).";
            free_idx.swap(zoned);
            if (free_idx.empty())
                logError() << "run_minimization : ZONE VIDE (aucune cellule value<0.5) !"
                           << " Verifier la carte de zone.";
        }
    }

    PorosityProblem prob(this, checkpoint_stride, alpha_reg,
                         phi_background.data(), free_idx, zonation);
    if (zonation) {
        logInfo() << "ZONATION active : " << free_idx.size() << " cellules libres groupees en "
                  << prob.n_zones << " zone(s) par valeur d'ebauche identique (chaque zone = 1 parametre).";
    }

    // --- Carte VERITE (porosity_ref_map) pour le suivi de l'ecart parametrique ---
    // sse_phi_ref = sum_i (phi_i - phi_ref_i)^2 ecrit dans descent_log.csv. Charge une
    // seule fois ; si absente/illisible, la colonne reste a -1 (non disponible).
    if (!porosity_ref_map_path.empty()) {
        std::vector<double> pref(cell_nb);
        const std::string ref_path = m_io->absolutePath(porosity_ref_map_path);
        if (load_phi_csv(ref_path, pref.data())) {
            prob.phi_ref = std::move(pref);
            logInfo() << "Suivi ecart parametrique sse_phi_ref vs " << ref_path;
        } else {
            logWarn() << "porosity_ref_map illisible (" << ref_path
                         << ") : colonne sse_phi_ref desactivee.";
        }
    } else {
        logInfo() << "porosity_ref_map non renseignee : colonne sse_phi_ref = -1.";
    }

    // --- Dossier de suivi de descente : <dossier>/descent/ (CSV + cartes phi/iteree) ---
    // Recree proprement a chaque run.
    {
        const std::string ddir = m_io->absolutePath("descent");
        std::error_code ec;
        std::filesystem::remove_all(ddir, ec);                 // repart propre
        std::filesystem::create_directories(ddir + "/phi_iter", ec);
        prob.descent_dir = ddir;
        logInfo() << "Suivi de descente -> " << ddir
                  << " (descent_log.csv + phi_iter/phi_iter_XXXX.csv)";
    }

    // Mode normalize : les bornes lo/hi sont calculees sur les cellules LIBRES.
    // On donne le masque au modele pour update_phi_norm_bounds().
    phi_norm_free_idx = free_idx;

    // Variable d'optimisation : theta = m (la porosite brute, non bornee).
    // On seede m par l'INVERSE du transform actif (swSpE::m_from_phi) pour que
    // phi_correct(theta) ~= phi_start :
    //   none/clamp : m = phi_start (identite sur la zone affine) ;
    //   tanh       : m = logit(phi_start), clampe pour eviter logit(1)=+inf ;
    //   normalize  : m = phi_start (les bornes lo/hi sont fixees ensuite par
    //                update_phi_norm_bounds, qui rescale de maniere coherente).
    // Sans cet inverse, en tanh la carte de depart serait silencieusement
    // phi=sigma(phi_start) (ex: phi_start=1 -> sigma(1)=0.731), pas phi_start.
    std::vector<double> theta((size_t)prob.size());
    if (zonation) {
        for (size_t z = 0; z < prob.n_zones; z++) theta[z] = m_from_phi(prob.zone_value[z]);
    } else {
        for (size_t j = 0; j < free_idx.size(); j++) theta[j] = m_from_phi(phi_start[free_idx[j]]);
    }

    const char* mode_name =
        (phi_transform == PhiTransform::Clamp)     ? "clamp (phi=clamp(m,0,1.05))" :
        (phi_transform == PhiTransform::Tanh)      ? "tanh (phi=1/(1+exp(-m)))" :
        (phi_transform == PhiTransform::Normalize) ? "normalize (phi=(m-lo)/(hi-lo), lo/hi=min/max cellules libres)" :
                                                     "none (phi=m brut)";
    const bool use_adam   = (minimize_algo == "adam");
    const bool use_gd     = (minimize_algo == "gd");
    const bool use_lbfgsb = (minimize_algo == "lbfgsb");
    logInfo() << "BEGIN MINIMIZATION ("
              << (use_gd ? "GD projetee + Armijo"
                         : (use_adam ? "Adam (optim)"
                                     : (use_lbfgsb ? "L-BFGS-B box (LBFGSpp)" : "PSQN BFGS dense")))
              << ", phi_transform + masque)";
    logInfo() << "  parametres : " << free_idx.size() << " cellules libres / " << cell_nb
              << ", max_it=" << max_it << ", checkpoint_stride=" << checkpoint_stride
              << ", phi_transform=" << mode_name
              << ", minimize_algo=" << minimize_algo
              << ", alpha_reg=" << alpha_reg;

    double J_final  = 0.0;   // cout J atteint a la fin (data + ridge)
    long   n_grad_done = 0;  // nb d'evaluations de gradient (forward+adjoint)

    if (use_gd)
    {
        // ===================== GD PROJETEE + ARMIJO (maison) =====================
        // POURQUOI cet algo (diagnostic 2026-06-29) : l'Adam "maison" normalise le
        // gradient par coordonnee (update_i = a * mhat_i/(sqrt(vhat_i)+eps) ~= a*sign(g_i)),
        // ce qui DETRUIT la magnitude du gradient. Or pour reconstruire un barrage la
        // descente est CONCENTREE sur quelques cellules (|g|~1e3 au barrage, ~0 ailleurs) :
        // etaler un pas egal sur les 5684 cellules libres donne une direction quasi inutile
        // -> J oscille/diverge quel que soit le pas (constate a step=1.0/0.02/0.005/0.001).
        // BFGS (PSQN) part lui d'un pas 1.0 * gradient(~2e3) -> overshoot massif (+544% J).
        //
        // Le gradient corrige (test de Taylor ~1, cf [[taylor-phiLR-hard-min]]) fournit une
        // VRAIE direction de descente -grad. On l'utilise telle quelle (magnitude respectee)
        // avec une recherche lineaire d'Armijo qui (a) auto-echelonne le pas malgre le
        // mauvais conditionnement, (b) GARANTIT J decroissant, (c) rejette tout pas NaN/inf
        // (typiquement quand phi->0 raidit les facteurs 1/phi du forward). La projection
        // clippe m dans [0,1] (modes none/clamp), comme l'Adam projete.
        //
        // Cout : 1 grad (build+adjoint+replay) / iteree + quelques prob.func (replay seul,
        // sur les MEMES records figes que le gradient -> coherent avec la pente d'Armijo).
        const bool   project = (phi_transform == PhiTransform::None ||
                                phi_transform == PhiTransform::Clamp);
        const double c1      = 1e-4;   // constante de decroissance suffisante (Armijo)
        const double shrink  = 0.5;    // facteur de backtracking
        const double grow    = 4.0;    // on retente un pas plus grand a chaque iteree
        const int    max_bt  = 40;     // nb max de backtracks par iteree
        const size_t nf      = (size_t)prob.size();   // = n_zones si zonation

        logInfo() << "  GD projetee + Armijo : pas initial=" << adam_step_size
                  << ", c1=" << c1 << ", shrink=" << shrink << ", iter_max=" << max_it
                  << ", projection m in [0,1]=" << (project ? "oui" : "non")
                  << " (phi_transform=" << mode_name << ")";

        std::vector<double> g(nf, 0.0), theta_trial(nf, 0.0);
        std::vector<double> best_theta(theta);
        double best_J  = std::numeric_limits<double>::max();
        long   best_it = 0, n_grad_local = 0;
        double t = (adam_step_size > 0.0) ? adam_step_size : 1e-6;  // graine de pas (auto-echelonnee)

        for (size_t it = 1; it <= max_it; ++it) {
            // Gradient + vrai cout J au point courant (prob.grad reconstruit les records a theta).
            double J = prob.grad(theta.data(), g.data());
            ++n_grad_local;
            if (J < best_J) { best_J = J; best_theta = theta; best_it = (long)it; }

            double gnorm2 = 0.0;
            for (size_t j = 0; j < nf; ++j) gnorm2 += g[j] * g[j];
            if (!(gnorm2 > 0.0)) { std::cout << "[GD] gradient nul/NaN -> arret (iter " << it << ").\n"; break; }

            // --- Recherche lineaire d'Armijo PROJETEE sur le cout de substitution (records figes) ---
            double t_try = t * grow;
            bool   accepted = false;
            double J_new = J;
            for (int b = 0; b < max_bt; ++b) {
                double pg = 0.0;  // pente realisee dans la direction -grad (tient compte de la projection)
                for (size_t j = 0; j < nf; ++j) {
                    double v = theta[j] - t_try * g[j];
                    if (project) { if (v < 0.0) v = 0.0; else if (v > 1.0) v = 1.0; }
                    theta_trial[j] = v;
                    pg += g[j] * (theta[j] - v);   // = sum g*(theta-theta_trial) >= 0
                }
                J_new = prob.func(theta_trial.data());
                if (std::isfinite(J_new) && J_new <= J - c1 * pg && pg > 0.0) { accepted = true; break; }
                t_try *= shrink;
            }

            if (!accepted) {
                std::cout << "[GD] pas de decroissance trouvee a l'iteree " << it
                          << " (pas mini=" << std::scientific << std::setprecision(3) << t_try
                          << ") -> arret. Meilleur itere #" << best_it << " (J=" << best_J << ").\n";
                break;
            }
            theta = theta_trial;   // pas accepte (J de substitution decroit)
            t     = t_try;         // on repart du pas accepte (re-amplifie de 'grow' a l'iteree suivante)
            std::cout << "[GD] iter " << it << "  J=" << std::scientific << std::setprecision(6) << J
                      << " -> J_surrogate=" << J_new << "  pas=" << t_try
                      << "  ||g||^2=" << gnorm2 << std::endl;
        }
        theta       = best_theta;
        J_final     = best_J;
        n_grad_done = n_grad_local;
        std::cout << "[GD] termine : meilleur J=" << std::scientific << std::setprecision(10)
                  << best_J << " (iteree #" << best_it << ")  n_grad=" << n_grad_local << std::endl;
    }
    else if (use_lbfgsb)
    {
        // ===================== L-BFGS-B box (LBFGSpp) =====================
        // L-BFGS-B impose nativement les bornes de boite sur la variable d'optimisation
        // theta (= m sur les cellules libres). Comme la box doit s'appliquer a PHI, on
        // travaille en transform "none" (phi = m) : la box [phi_lower_bound, phi_upper_bound]
        // sur theta EST alors la box sur phi. Avec un autre transform, theta != phi et la
        // box ne porterait pas sur phi -> on avertit.
        if (phi_transform != PhiTransform::None) {
            logWarn() << "  L-BFGS-B : phi_transform=" << mode_name
                      << " != none -> la box [" << phi_lower_bound << ", " << phi_upper_bound
                      << "] porte sur m (var. brute), PAS directement sur phi. "
                      << "Recommande : phi_transform none.";
        }

        using LbVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
        const size_t nf = (size_t)prob.size();   // = n_zones si zonation, sinon nb cellules libres

        // Adaptateur objectif pour LBFGSpp : operator()(x, grad) -> f, en deleguant a
        // prob.grad (qui calcule J = J_data + ridge + L1, le gradient sur les libres, ET
        // ecrit la ligne descent_log + la carte phi a chaque evaluation).
        struct LbfgsbObjective {
            PorosityProblem& prob;
            double operator()(const LbVector& x, LbVector& grad) {
                return prob.grad(x.data(), grad.data());
            }
        } objfun{prob};

        LBFGSpp::LBFGSBParam<double> lbparam;
        lbparam.m              = lbfgsb_m;
        lbparam.epsilon        = lbfgsb_epsilon;
        lbparam.epsilon_rel    = lbfgsb_epsilon;
        lbparam.max_iterations = (int)max_it;
        lbparam.max_linesearch = lbfgsb_max_linesearch;

        logInfo() << "  L-BFGS-B parameters : box phi in [" << phi_lower_bound << ", "
                  << phi_upper_bound << "], m=" << lbparam.m
                  << ", epsilon=" << lbparam.epsilon
                  << ", max_linesearch=" << lbparam.max_linesearch
                  << ", max_iterations=" << lbparam.max_iterations;

        LBFGSpp::LBFGSBSolver<double> solver(lbparam);

        LbVector x(nf), lb(nf), ub(nf);
        for (size_t j = 0; j < nf; ++j) {
            x[j]  = theta[j];
            lb[j] = phi_lower_bound;
            ub[j] = phi_upper_bound;
        }

        // FIX gel L-BFGS-B : une evaluation "sonde" au point de depart fixe le segment
        // d'explosion de l'adjoint ; on FIGE ensuite la fenetre de troncature a ce segment
        // pour tout le solve (forced_first_seg_), afin que chaque evaluation de line-search
        // calcule J sur la MEME fenetre. Sans cela, le point d'explosion (fonction de theta)
        // fait varier la fenetre de J entre 2 essais -> Wolfe echoue -> l'optimiseur gele.
        {
            LbVector g_probe(nf);
            objfun(x, g_probe);                        // pose adjoint_truncation_seg_
            forced_first_seg_ = std::max(0, adjoint_truncation_seg_ + 1);
            logInfo() << "  L-BFGS-B : fenetre de troncature figee au segment "
                      << forced_first_seg_ << " (coherence line-search).";
        }

        double fx = std::numeric_limits<double>::quiet_NaN();
        int niter = 0;
        try {
            niter = solver.minimize(objfun, x, fx, lb, ub);
        } catch (const std::exception& e) {
            logWarn() << "  L-BFGS-B interrompu (" << e.what()
                      << ") : on garde le dernier itere atteint.";
        }
        forced_first_seg_ = -1;   // restaure la troncature dynamique (autres modes)

        for (size_t j = 0; j < nf; ++j) theta[j] = x[j];
        J_final     = fx;
        n_grad_done = (long)prob.grad_calls;   // 1 eval (forward+adjoint) par appel objfun
        std::cout << "[L-BFGS-B] termine : J_final=" << std::scientific << std::setprecision(10)
                  << fx << "  iterations=" << niter
                  << "  ||proj grad||_inf=" << solver.final_grad_norm()
                  << "  n_eval(grad)=" << prob.grad_calls << std::endl;
    }
    else if (!use_adam)
    {
        // ===================== BFGS (PSQN) =====================
        // BFGS dense : memoire de travail = 7*n + n*(n+1)/2 doubles (alloc interne).
        PSQN::optim_info res = PSQN::bfgs<MinimizeReporter>(
            prob, theta.data(),
            /*rel_eps*/ 1e-8,
            /*max_it */ (PSQN::psqn_uint)max_it,
            /*c1     */ 1e-4,
            /*c2     */ 0.9,
            /*trace  */ 1);
        J_final     = res.value;
        n_grad_done = (long)res.n_grad;
        std::cout << "[BFGS] code=" << (int)res.info
                  << "  J_final=" << std::scientific << std::setprecision(10) << res.value
                  << "  n_eval=" << res.n_eval << "  n_grad=" << res.n_grad << std::endl;
    }
    else
    {
        // ===================== Adam PROJETE (maison) =====================
        // Pourquoi pas optim::gd : optim impose les bornes par un TRANSFORM interne
        // (box_objfn -> sigmoide), ce qui REINTRODUIT la saturation dphi/dm->0 qu'on
        // veut justement eviter. On ecrit donc Adam (method 6) a la main + PROJECTION :
        // apres chaque update, m est clippe dans [0,1] quand le transform rend m == phi
        // (modes none/clamp). Avantage vs clamp-dans-la-tape : une cellule SUR la borne
        // dont le gradient pointe vers l'interieur REPART (le clip ne bloque que la
        // composante sortante) -> pas de gel aveugle, et phi=0 atteignable exactement.
        //
        // Adam (sign-like par coordonnee) : update_i = a * mhat_i / (sqrt(vhat_i)+eps).
        // Toutes les cellules (dont le barrage) avancent de ~a, independamment de la
        // magnitude du gradient. Chaque iteree = 1 forward + 1 adjoint pave (prob.grad).
        const double a   = adam_step_size;
        const double b1  = 0.9, b2 = 0.999, eps = 1e-8;
        // Projection de m dans [0,1] seulement quand m EST la porosite (none/clamp).
        const bool   project = (phi_transform == PhiTransform::None ||
                                phi_transform == PhiTransform::Clamp);

        logInfo() << "  Adam parameters : step=" << a << ", beta_1=" << b1
                  << ", beta_2=" << b2 << ", iter_max=" << max_it
                  << ", projection m in [0,1]=" << (project ? "oui" : "non")
                  << " (phi_transform=" << mode_name << ")";

        // Garde anti-explosion adjointe : , 
        // On surveille ||g||_inf : un saut > EXPLOSION_FACT * (iteree precedente) signe
        // une explosion -> on ARRETE (sans appliquer le pas pourri) et on garde le
        // MEILLEUR itere (J min). On ecrit donc toujours une carte saine, jamais la divergee.
        const double EXPLOSION_FACT = 50.0;
        const size_t nf = (size_t)prob.size();   // = n_zones si zonation
        std::vector<double> g(nf, 0.0), mom(nf, 0.0), vel(nf, 0.0);
        std::vector<double> best_theta(theta);
        double best_J = std::numeric_limits<double>::max();
        long   best_it = 0;
        double g_inf_prev = -1.0;
        long n_grad_local = 0;
        for (size_t it = 1; it <= max_it; ++it) {
            // Gradient au point courant (prob.grad ecrit aussi la ligne de suivi + la carte phi).
            J_final = prob.grad(theta.data(), g.data());
            ++n_grad_local;

            // Suivi du meilleur itere (le cout J vient du forward, fiable meme si l'adjoint
            // se degrade ; on garde la carte qui a donne le plus petit J).
            if (J_final < best_J) { best_J = J_final; best_theta = theta; best_it = (long)it; }

            // Detection d'explosion : saut soudain de ||g||_inf. STOP UNIQUEMENT si
            // adam_safe=true (defaut false). Sinon on signale mais on CONTINUE : les pics
            // transitoires se resorbent souvent et la descente reprend (le meilleur-itere
            // protege de toute facon la carte finale).
            double g_inf = 0.0;
            for (size_t j = 0; j < nf; ++j) g_inf = std::max(g_inf, std::fabs(g[j]));
            if (g_inf_prev > 0.0 && g_inf > EXPLOSION_FACT * g_inf_prev) {
                std::cout << "[Adam] " << (adam_safe ? "*** EXPLOSION ADJOINTE -> ARRET (adam_safe)"
                                                     : "(pic ||g|| detecte, on continue : adam_safe=false)")
                          << " a l'iteree " << it
                          << " : ||g||_inf=" << std::scientific << std::setprecision(3) << g_inf
                          << " >> " << EXPLOSION_FACT << "x precedent (" << g_inf_prev << ")."
                          << "  meilleur itere #" << best_it << " (J=" << best_J << ").\n";
                if (adam_safe) break;
            }
            g_inf_prev = g_inf;

            // Update Adam projete.
            const double bc1 = 1.0 - std::pow(b1, (double)it);   // bias-correction
            const double bc2 = 1.0 - std::pow(b2, (double)it);
            for (size_t j = 0; j < nf; ++j) {
                mom[j] = b1 * mom[j] + (1.0 - b1) * g[j];
                vel[j] = b2 * vel[j] + (1.0 - b2) * g[j] * g[j];
                const double mhat = mom[j] / bc1;
                const double vhat = vel[j] / bc2;
                theta[j] -= a * mhat / (std::sqrt(vhat) + eps);   // descente
                if (project) {                                    // clip dans [0,1]
                    if (theta[j] < 0.0) theta[j] = 0.0;
                    else if (theta[j] > 1.0) theta[j] = 1.0;
                }
            }
        }
        // On retourne le MEILLEUR itere (pas le dernier, qui peut avoir diverge).
        theta   = best_theta;
        J_final = best_J;
        n_grad_done = n_grad_local;
        std::cout << "[Adam] termine : meilleur J=" << std::scientific << std::setprecision(10)
                  << best_J << " (iteree #" << best_it << ")  n_grad=" << n_grad_local
                  << "  projection=" << (project ? "[0,1]" : "off") << std::endl;
    }

    logInfo() << "END MINIMIZATION";

    taylor_test_after_gradient = true;  // restore default

    // theta(=m) optimise -> phi physique : cellules fixees = phi_b, libres = phi_correct(m).
    // (meme transform que dans la tape -> coherence parfaite avec le cout optimise.)
    // Pour normalize, on refige lo/hi sur la carte FINALE (sa propre min/max) avant
    // de transformer, comme pour les iterees precedentes.
    // Expansion zone-aware : theta est par-cellule (defaut) OU par-zone (zonation). On
    // recupere la valeur m de chaque cellule libre via l'indirection de zone si active.
    auto theta_of_free = [&](size_t j) -> double {
        return zonation ? theta[prob.zone_of_free[j]] : theta[j];
    };
    std::vector<double> m_final(cell_nb);
    for (size_t i = 0; i < cell_nb; i++) m_final[i] = phi_background[i];
    for (size_t j = 0; j < free_idx.size(); j++) m_final[free_idx[j]] = theta_of_free(j);
    update_phi_norm_bounds(m_final.data());

    for (size_t i = 0; i < cell_nb; i++) phi_start[i] = phi_background[i];
    for (size_t j = 0; j < free_idx.size(); j++) {
        phi_start[free_idx[j]] = phi_correct(theta_of_free(j));
    }

    // Decomposition finale : J_total = J_data + (alpha/2)||phi - phi_b||^2
    double reg_norm2 = 0.0;
    for (size_t i = 0; i < cell_nb; i++) {
        double d = phi_start[i] - phi_background[i];
        reg_norm2 += d * d;
    }
    double J_data_final = J_final - 0.5 * alpha_reg * reg_norm2;
    std::cout << "[Minim] alpha=" << alpha_reg
              << "  J_data=" << std::scientific << std::setprecision(10) << J_data_final
              << "  ||phi-phi_b||^2=" << reg_norm2
              << "  J_total=" << J_final << std::endl;

    // Ecriture de la carte optimisee au format SP_cell_porosity_map.txt (Unif/Distrib).
    // -> directement re-utilisable comme carte de porosite (cp vers input/SP_cell_porosity_map.txt)
    //    et traçable par les outils mesh_plot. Extension .txt (le contenu n'est plus du CSV).
    const std::string out_map = m_io->absolutePath("optimized_porosity.txt");
    write_sp_porosity_map(out_map, phi_start, cell_nb);
    logInfo() << "Carte de porosite optimisee (format SP, rechargeable tel quel comme "
              << "SP_cell_porosity_map.txt) ecrite dans : " << out_map;

    // Ligne pour la courbe en L : <dossier>/alpha_sweep.csv (cree avec en-tete si absent)
    const std::string sweep_csv = m_io->absolutePath("alpha_sweep.csv");
    bool sweep_exists = std::filesystem::exists(sweep_csv);
    std::ofstream sfs(sweep_csv, std::ios::app);
    if (sfs.is_open()) {
        if (!sweep_exists) sfs << "alpha,J_data,reg_norm2,J_total,n_grad\n";
        sfs << std::setprecision(12) << alpha_reg << "," << J_data_final << ","
            << reg_norm2 << "," << J_final << "," << n_grad_done << "\n";
        sfs.close();
        logInfo() << "Ligne ajoutee a la courbe en L : " << sweep_csv;
    }

    return J_final;
}
