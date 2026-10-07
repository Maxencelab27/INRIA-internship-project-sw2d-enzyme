/**
 * Checkpointed adjoint computation for porosity calibration.
 *
 * Memory model
 * ============
 * Without checkpointing, Enzyme tapes every intermediate value for every
 * time step → O(T * mesh_size) memory → OOM on large cases.
 *
 * With checkpointing (stride S):
 *   - Tape per Enzyme call : O(S * mesh_size)   ← bounded, fits in RAM
 *   - Checkpoint storage   : O((T/S) * 3*N_cells) ← tiny (no edge arrays)
 *   - Recomputation cost   : (T/S) extra forward passes of S steps each
 *
 * For a mesh with 50 000 cells / 150 000 edges and 2 000 time steps:
 *   Without : ~120 GB tape  → crash
 *   With S=30: ~1.8 GB tape + ~150 MB checkpoints → runs fine
 *
 * Implementation strategy
 * =======================
 * 1. build_checkpoint_records()  – dry forward run; records dt, BC values,
 *    and observation triggers per step; saves primal state every S steps.
 *
 * 2. ckpt_segment_E()  – static free function that Enzyme differentiates.
 *    Takes prerecorded dt/BC as enzyme_const raw pointers (they are read-only
 *    data, Enzyme must not track adjoints through them).
 *    Runs S physics steps and returns the accumulated cost for those steps.
 *
 * 3. get_spatial_gradient_checkpointed()  – backward loop over segments.
 *    For k = K-1 … 0:
 *      a) Restore primal solver from checkpoint k.
 *      b) Seed adjoint_solver's state with d_h / d_qx / d_qy from segment k+1.
 *      c) Call __enzyme_autodiff on ckpt_segment_E.
 *      d) Harvest new d_h / d_qx / d_qy (adjoint of input state) for segment k-1.
 *      e) dJ_dPhi accumulates via the enzyme_dup phi_map / dJ_dPhi pair.
 */

#include "swSpE.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <filesystem>
#include "../../../common/swLog.h"

// ---------------------------------------------------------------------------
// Enzyme interface
// ---------------------------------------------------------------------------
extern int enzyme_dup;
extern int enzyme_const;
extern int enzyme_out;
extern int enzyme_dupnoneed;

extern "C" {
    void __enzyme_autodiff(void*, ...);
}

// ---------------------------------------------------------------------------
// Sonde mémoire optionnelle (activée par SW2D_MEM_PROBE=1). Sert à diagnostiquer
// la croissance RAM du forward _E en gen-ref-e :
//   - vmrss_kb()         : RSS réelle du process (pages réellement résidentes).
//   - heap_inuse_bytes() : octets HEAP réellement *en cours d'utilisation*
//     (mallinfo2.uordblks + hblkhd). Si ça monte régulièrement -> VRAIE fuite ;
//     si RSS monte mais inuse reste plat (et malloc_trim fait retomber la RSS)
//     -> simple fragmentation d'arène glibc, pas une fuite.
// Aucun coût quand la variable d'env n'est pas définie (tout est gardé).
// ---------------------------------------------------------------------------
#include <malloc.h>
namespace {
    inline long vmrss_kb() {
        std::ifstream f("/proc/self/status");
        std::string key; long val = -1;
        while (f >> key) {
            if (key == "VmRSS:") { f >> val; break; }
        }
        return val;
    }
    inline size_t heap_inuse_bytes() {
#if defined(__GLIBC__) && (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 33))
        struct mallinfo2 mi = mallinfo2();
        return (size_t)mi.uordblks + (size_t)mi.hblkhd;
#else
        struct mallinfo mi = mallinfo();   // attention: champs int -> peut déborder >2 Go
        return (size_t)(unsigned)mi.uordblks + (size_t)(unsigned)mi.hblkhd;
#endif
    }
}

// ===========================================================================
// 1.  build_checkpoint_records
//     Dry-run: full forward pass that records per-step data and checkpoint
//     states. Does NOT use Enzyme at all — just plain physics.
// ===========================================================================
void swSpE::build_checkpoint_records(double* phi_map, size_t stride)
{
    segs_record_.clear();

    // --- Initialise primal state ---
    this->reset_to_start();
    // Changement de variable : phi_map est la variable d'optimisation m ; phi = phi_correct(m)
    // (none/clamp/tanh/normalize, pilote par phi_transform). phi_ptr pointe sur le phi
    // TRANSFORME (phi.data()), pas sur m (phi_map). Les bornes de normalisation sont figees
    // en amont (update_phi_norm_bounds, dans get_spatial_gradient_checkpointed).
    for (size_t i = 0; i < cell_nb; i++) phi[i] = phi_correct(phi_map[i]);
    this->phi_ptr = phi.data();

    // FIX A : les facteurs d'arête (1/phi) sont precalcules a l'init avec la carte
    // d'entree. Comme phi change ici (assimilation), il FAUT les recalculer, sinon
    // le bilan applique des flux incoherents (ex: 1/phi_true=x50 au barrage) -> q
    // explose. NB : ce recalcul reste HORS de la tape Enzyme (computeEdgeFactors
    // n'est pas differencie) -> le gradient ne capte pas encore le terme 1/phi
    // (bug C, a traiter via un computeEdgeFactors_E). Mais le forward redevient
    // physiquement coherent pour le phi courant.
    this->computeEdgeFactors();

    // ---- Dump one-shot des cellules de bord + type BC (SW2D_DUMP_BND=1) ----
    // Permet de savoir si une cellule donnée (ex: cell 183 qui s'emballe) est sur un bord
    // et de quel type (q / fr / h / z), pour localiser un bug de BC vs un bug de schéma.
    if (std::getenv("SW2D_DUMP_BND")) {
        std::ofstream o("/tmp/E_bnd_cells.csv");
        o << "edge\tcell\tbc_type\tbc_val_t0\tx\ty\tzb\tK1\tK2\n";
        o << std::setprecision(8);
        for (size_t e = 0; e < edge_bnd_nb; e++) {
            size_t c = edge_bnd_cell_idx[e];
            std::string bt = edge_bnd_bc_pointer[e] ? edge_bnd_bc_pointer[e]->getBCType() : "?";
            double k1 = friction_coeff_direction_1.empty() ? -1.0 : friction_coeff_direction_1[c];
            double k2 = friction_coeff_direction_2.empty() ? -1.0 : friction_coeff_direction_2[c];
            o << e << '\t' << c << '\t' << bt << '\t'
              << (e < bc_val_interpolated.size() ? bc_val_interpolated[e] : 0.0) << '\t'
              << cell_centre_x[c] << '\t' << cell_centre_y[c] << '\t' << cell_z_b[c]
              << '\t' << k1 << '\t' << k2 << '\n';
        }
        // friction des cellules suspectes (toutes) -> repérer K=0 non initialisés
        std::ofstream of("/tmp/E_friction.csv");
        of << "cell\tx\ty\tK1\tK2\n" << std::setprecision(8);
        for (size_t c = 0; c < cell_nb; c++) {
            double k1 = friction_coeff_direction_1.empty() ? -1.0 : friction_coeff_direction_1[c];
            double k2 = friction_coeff_direction_2.empty() ? -1.0 : friction_coeff_direction_2[c];
            of << c << '\t' << cell_centre_x[c] << '\t' << cell_centre_y[c] << '\t' << k1 << '\t' << k2 << '\n';
        }
        std::cout << "[DUMP_BND] " << edge_bnd_nb << " bnd cells -> /tmp/E_bnd_cells.csv ; friction -> /tmp/E_friction.csv\n";
    }

    size_t obs_idx      = 0;
    size_t step_in_seg  = 0;

    SegPhysData cur;
    cur.t_start  = this->t;
    cur.n_steps  = 0;
    cur.h_ckpt   = water_depth;
    cur.qx_ckpt  = unit_discharge_x;
    cur.qy_ckpt  = unit_discharge_y;

    // ---- Sonde mémoire (SW2D_MEM_PROBE=1) ----
    const bool   mem_probe    = (std::getenv("SW2D_MEM_PROBE") != nullptr);
    const size_t probe_stride = 50;     // résumé RSS/inuse tous les 50 pas
    size_t       step_global  = 0;
    size_t       inuse_ref    = mem_probe ? heap_inuse_bytes() : 0; // inuse au pas 0
    if (mem_probe) {
        std::cout << "[MEMPROBE] init: rss=" << vmrss_kb()/1024.0 << "MB inuse="
                  << inuse_ref/1048576.0 << "MB  cells=" << cell_nb
                  << " edge_int=" << edge_int_nb << " edge_bnd=" << edge_bnd_nb
                  << "  (segs_record_ ~ " << (3*cell_nb*8 + 30*edge_bnd_nb*8)/1024
                  << "KB/segment)\n";
    }

    // ---- Dump de champ aux temps cibles (SW2D_DUMP_TIMES="400,450,460,463,465") ----
    // Écrit h/qx/qy par cellule du forward _E pour comparaison avec le solveur classique
    // (mêmes colonnes que ResultHydro_*.txt). Sert à localiser où le _E diverge à t≈466 s.
    std::vector<double> dump_times;
    if (const char* dt_env = std::getenv("SW2D_DUMP_TIMES")) {
        std::stringstream ss(dt_env); std::string tok;
        while (std::getline(ss, tok, ',')) { try { dump_times.push_back(std::stod(tok)); } catch (...) {} }
        std::sort(dump_times.begin(), dump_times.end());
    }
    size_t next_dump = 0;

    // ---- Audit du clamp h<0 (SW2D_CLAMP_AUDIT=1) ----
    // Le bilan _E peut produire des profondeurs negatives (overshoot de flux) ; le clamp
    // h<0 -> h=0 RECREE alors de la masse. On cumule ce volume (m^3, via cell_area) et le
    // nombre d'evenements pour voir si l'exces de masse vient de la (vs un desequilibre de flux).
    const bool clamp_audit = (std::getenv("SW2D_CLAMP_AUDIT") != nullptr);
    double clamp_mass_created = 0.0;   // volume total recree par le clamp (m^3)
    double clamp_mass_neg_max = 0.0;   // pire h negatif rencontre (m)
    size_t clamp_events = 0;           // nb total de cellules clampees
    double bnd_inflow_total = 0.0;     // masse cumulee injectee par TOUTES les aretes de bord (m^3)
    double bnd_inflow_q_total = 0.0;   // idem mais seulement aretes type "q" (entree)
    size_t q_crit_branch = 0, q_raw_branch = 0; // nb d'aretes "q" en regime critique vs bc_val brut
    auto current_volume = [&]() {
        double v = 0.0; for (size_t i = 0; i < cell_nb; i++) v += water_depth[i] * cell_area[i];
        return v;
    };

    auto dump_field = [&](double tgt) {
        std::ostringstream fn; fn << "/tmp/E_dump_t" << (int)std::lround(tgt) << ".csv";
        std::ofstream o(fn.str());
        o << "cell\tx\ty\tzb\th\tqx\tqy\tarea\n";
        o << std::setprecision(10);
        for (size_t i = 0; i < cell_nb; i++)
            o << i << '\t' << cell_centre_x[i] << '\t' << cell_centre_y[i] << '\t'
              << cell_z_b[i] << '\t' << water_depth[i] << '\t'
              << unit_discharge_x[i] << '\t' << unit_discharge_y[i] << '\t'
              << cell_area[i] << '\n';
        std::cout << "[DUMP_E] t=" << t << " (cible " << tgt << ") -> " << fn.str() << "\n";
    };

    while (t < t_max)
    {
        // Échantillon par-phase sur quelques pas ciblés (localise la fonction qui alloue)
        const bool mem_sample = mem_probe && (step_global == 100 || step_global == 500
                                              || step_global == 2000);
        size_t mph = mem_sample ? heap_inuse_bytes() : 0;
        auto phase = [&](const char* name) {
            if (!mem_sample) return;
            size_t now = heap_inuse_bytes();
            std::cout << "[MEMPHASE step=" << step_global << "] " << name << " d_inuse="
                      << (long long)(now - mph) << " B\n";
            mph = now;
        };

        // ---- Precompute BC at current t (stores into bc_val_interpolated) ----
        computeBC_E(t);
        phase("computeBC_E");
        // Copy BC snapshot BEFORE the step (this is what hypComputeFluxes_E reads)
        std::vector<double> bc_snap(bc_val_interpolated.begin(),
                                     bc_val_interpolated.end());

        // ---- Decide target dt (clamp to next observation) ----
        double dt_target = dt_max;
        bool   at_obs    = false;
        if (obs_idx < current_n_obs) {
            double dt_to_obs = current_obs_times[obs_idx] - t;
            if (dt_to_obs - 10e-6 > 0.0 && dt_to_obs < dt_target)
                dt_target = dt_to_obs;
        }

        // ---- MUSCL LR-states + CFL (same order as hypCompute_E) ----
        // hypCompute_E does: LRStates → TimeStep → Fluxes → Balance.
        // We replicate that here so lambda arrays are populated before TimeStep reads them.
        hypComputeLRStates_E();
        phase("hypComputeLRStates_E");
        double dt_cfl = hypComputeTimeStep_E();
        phase("hypComputeTimeStep_E");

        double cur_dt = (dt_cfl < dt_target) ? dt_cfl : dt_target;
        if ((t + cur_dt) > t_max) cur_dt = t_max - t;
        if (cur_dt < 1e-6) cur_dt = 1e-6;
        this->dt = cur_dt;

        // ---- Run one physics step ----
        hypComputeFluxes_E();
        phase("hypComputeFluxes_E");
        hypComputeBalance_E(cur_dt);
        phase("hypComputeBalance_E");
        for (size_t n = 0; n < cell_nb; ++n) {
            if (water_depth[n] < 0.) {
                if (clamp_audit) {
                    clamp_mass_created += (-water_depth[n]) * cell_area[n];
                    if (-water_depth[n] > clamp_mass_neg_max) clamp_mass_neg_max = -water_depth[n];
                    clamp_events++;
                }
                water_depth[n] = 0.;
                unit_discharge_x[n] = 0.;
                unit_discharge_y[n] = 0.;
            }
        }
        if (clamp_audit) {
            // Masse injectée par les arêtes de bord ce pas (= ce que hypComputeBalance_E a
            // ajouté côté bord). Sépare aussi le débit "q" et compte le régime critique.
            for (size_t e = 0; e < edge_bnd_nb; e++) {
                double m = edge_bnd_factor_1_F1[e] * F1_bnd[e] * cur_dt
                         * cell_area[edge_bnd_cell_idx[e]];
                bnd_inflow_total += m;
                std::string bt = edge_bnd_bc_pointer[e] ? edge_bnd_bc_pointer[e]->getBCType() : "";
                if (bt == "q") {
                    bnd_inflow_q_total += m;
                    double bcv = (e < bc_val_interpolated.size()) ? bc_val_interpolated[e] : 0.0;
                    // régime critique si F1_bnd s'écarte nettement du bc_val brut
                    if (bcv > small_real && std::abs(F1_bnd[e] - bcv) > 0.05 * bcv) q_crit_branch++;
                    else if (bcv > small_real) q_raw_branch++;
                }
            }
        }
        // Friction : le classique ne l'applique QUE si model_friction != "none"
        // (execute(): `if (model_friction != "none") sourceFriction(...)`). Avec le défaut
        // model_friction="none", AUCUNE friction côté classique. Appeler
        // sourceFrictionStrickler_E inconditionnellement injectait une friction fantôme
        // (K=20 par défaut car friction_coeff non initialisé) -> ralentissait la propagation
        // -> sous-mouillage vs forward classique de référence.
        if (model_friction != "none")
            sourceFrictionStrickler_E(cur_dt);
        phase("sourceFrictionStrickler_E");
        // Head loss singulier APRES la friction (meme ordre que execute() classique).
        // Doit etre IDENTIQUE a ckpt_segment_E (coherence dry-run/tape) sinon les
        // checkpoints sont hors trajectoire -> gradient biaise. Son absence rendait le
        // forward _E non dissipatif -> adjoint explosif (cas Toce, Taylor casse).
        if (model_singular_head_loss != "none")
            sourceSingularHeadLoss_E(cur_dt);
        t += cur_dt;

        // ---- Dump de champ aux temps cibles (premier pas qui dépasse la cible) ----
        while (next_dump < dump_times.size() && t >= dump_times[next_dump]) {
            dump_field(dump_times[next_dump]);
            // Carte native ResultHydro_p*.txt (même format que le solveur classique) du
            // forward _E à l'instant cible, écrite dans <folder>/output/.
            storeResultMaps(map_storage_mode, t);
            if (clamp_audit) {
                double vol = current_volume();
                std::cout << "[CLAMP_AUDIT] t=" << t << " step=" << step_global
                          << "  volume=" << vol << " m3"
                          << "  masse_recreee_clamp(cumul)=" << clamp_mass_created << " m3"
                          << " (" << (vol > 0 ? 100.0*clamp_mass_created/vol : 0.0) << "% du volume)"
                          << "  events=" << clamp_events
                          << "  pire_h_neg=" << clamp_mass_neg_max << " m"
                          << "  | inflow_bord=" << bnd_inflow_total << " m3"
                          << " (dont q=" << bnd_inflow_q_total << " m3)"
                          << "  q_critique/brut=" << q_crit_branch << "/" << q_raw_branch
                          << std::endl;
            }
            next_dump++;
        }

        // [MASSDBG] diagnostic temporaire : volume d'eau total du forward _E (SW2D_MASS_DEBUG=1)
        {
            static const bool mass_dbg = (std::getenv("SW2D_MASS_DEBUG") != nullptr);
            if (mass_dbg && (t < 30.0 || std::fmod(t, 500.0) < cur_dt)) {
                double vol = 0.0, sumh = 0.0; int wet = 0;
                for (size_t i = 0; i < cell_nb; i++) {
                    double h = water_depth[i];
                    vol += h * cell_area[i]; sumh += h;
                    if (h > 0.01) wet++;
                }
                std::cout << "[MASSDBG E] t=" << t << " vol=" << vol
                          << " sumh=" << sumh << " wet=" << wet << "/" << cell_nb << std::endl;
            }
        }

        // ---- Check observation trigger ----
        // Mécanisme en deux temps :
        //
        // (a) SKIP : avancer obs_idx si l'observation a été dépassée sans déclencher.
        //     Cela évite qu'une obs à t=0 (dt_to_obs=0, jamais clampée) bloque
        //     toutes les observations suivantes.
        while (obs_idx < current_n_obs && t > current_obs_times[obs_idx] + 1e-6) {
            std::cout << "[Checkpoint] WARNING obs " << obs_idx
                      << " (t=" << current_obs_times[obs_idx] << "s) dépassée"
                      << " (t_sim=" << t << "s), saut automatique.\n";
            obs_idx++;
        }
        // (b) TRIGGER : le pas a atterri exactement sur le temps d'observation.
        //     On compare le TEMPS SIMULE t a l'obs, et NON cur_dt a dt_target.
        //     BUG corrige : dt_target vaut dt_max tant que l'obs est loin (clamp
        //     ligne ~96 seulement si dt_to_obs < dt_max). Des que la CFL cesse de
        //     limiter (dt_cfl >= dt_max), cur_dt == dt_max == dt_target et l'ancienne
        //     condition |cur_dt - dt_target| < 1e-11 declenchait TOUTES les obs en
        //     rafale au mauvais instant (ex: t~3900s au lieu de 13800-14400s).
        int triggered_obs = -1;
        if (obs_idx < current_n_obs &&
            std::abs(t - current_obs_times[obs_idx]) < 1e-6)
        {
            triggered_obs = (int)obs_idx;
            obs_idx++;
        }

        // ---- gen-ref (_E forward) : dump the reference state at each obs trigger ----
        // gen_ref_mode is set ONLY by generate_references_E (off during assimilation, so
        // this is a no-op there). We dump the CURRENT _E state (water_depth/qx/qy) at the
        // exact obs time. Assimilation later replays this SAME _E forward (same phi, same
        // dt sequence, same obs grid) -> water_depth matches the dumped reference -> J = 0.
        if ((gen_ref_mode || forward_dump_mode) && triggered_obs >= 0) {
            storeReferenceMap(t);   // ecrit dans state_dump_dir (reference/ ou output/)
        }

        // ---- Store step record ----
        cur.dt_seq.push_back(cur_dt);
        // Append BC snapshot for this step (flat: step * edge_bnd_nb + j)
        cur.bc_flat.insert(cur.bc_flat.end(), bc_snap.begin(), bc_snap.end());
        cur.obs_trig.push_back(triggered_obs);
        cur.n_steps++;
        step_in_seg++;

        // ---- Checkpoint boundary ----
        bool segment_done = (step_in_seg == stride) || (t >= t_max);
        if (segment_done) {
            segs_record_.push_back(cur); // copy-store complete segment

            step_in_seg = 0;
            if (t < t_max) {
                cur.t_start  = t;
                cur.n_steps  = 0;
                cur.h_ckpt   = water_depth;
                cur.qx_ckpt  = unit_discharge_x;
                cur.qy_ckpt  = unit_discharge_y;
                cur.dt_seq.clear();
                cur.bc_flat.clear();
                cur.obs_trig.clear();
            }
        }

        // ---- Sonde mémoire : résumé périodique + test fragmentation (malloc_trim) ----
        if (mem_probe && (step_global % probe_stride == 0)) {
            size_t inuse = heap_inuse_bytes();
            long   rss_before = vmrss_kb();
            // estimation de la part imputable à segs_record_ (le reste = "autre")
            size_t segs_bytes = 0;
            for (const auto& s : segs_record_)
                segs_bytes += s.h_ckpt.capacity()*8 + s.qx_ckpt.capacity()*8
                            + s.qy_ckpt.capacity()*8 + s.bc_flat.capacity()*8
                            + s.dt_seq.capacity()*8  + s.obs_trig.capacity()*4;
            malloc_trim(0);                      // rend au noyau ce qui est libre
            long rss_after = vmrss_kb();
            std::cout << "[MEMPROBE] step=" << step_global << " t=" << t
                      << "  rss=" << rss_before/1024.0 << "MB"
                      << "  inuse=" << inuse/1048576.0 << "MB"
                      << "  d_inuse_vs_step0=" << (long long)(inuse - inuse_ref)/1048576.0 << "MB"
                      << "  segs=" << segs_record_.size() << "(" << segs_bytes/1048576.0 << "MB)"
                      << "  autre=" << (long long)(inuse - inuse_ref - (long long)segs_bytes)/1048576.0 << "MB"
                      << "  | rss_apres_trim=" << rss_after/1024.0 << "MB"
                      << "  (delta_trim=" << (rss_before-rss_after)/1024.0 << "MB => fragmentation si >0)\n";
        }
        step_global++;
    }

    // --- Bilan forward : vérifie que toutes les observations ont bien été déclenchées ---
    int total_obs_triggered = 0;
    for (const auto& s : segs_record_)
        for (int oi : s.obs_trig)
            if (oi >= 0) total_obs_triggered++;

    if (!minim_quiet_logs) {
    std::cout << "[Checkpoint] " << segs_record_.size()
              << " segments recorded (stride=" << stride << ", total steps≈"
              << (segs_record_.empty() ? 0 :
                  segs_record_.back().n_steps + (segs_record_.size()-1)*stride)
              << ")\n";
    std::cout << "[Checkpoint] Observations déclenchées dans le forward run : "
              << total_obs_triggered << " / " << current_n_obs << "\n";
    }
    if (total_obs_triggered != (int)current_n_obs) {
        std::cout << "[Checkpoint] ATTENTION : "
                  << (current_n_obs - total_obs_triggered)
                  << " observation(s) non déclenchée(s) ! "
                  << "Vérifiez les temps dans observations_assimilation.dat "
                  << "(doivent être en secondes, < t_max=" << t_max << "s)\n";
        // Afficher les temps attendus vs t_max
        for (size_t i = 0; i < current_n_obs; i++) {
            std::cout << "  obs[" << i << "] t=" << current_obs_times[i] << "s"
                      << (current_obs_times[i] > t_max ? "  <- HORS FENETRE !" : "")
                      << "\n";
        }
    }
}


// ===========================================================================
// 1bis. generate_references_E
//     gen-ref via the _E forward, so J(phi_ref) == 0 by construction.
//
//     The classic gen-ref writes references from execute() (the NON-Enzyme engine),
//     while assimilation evaluates J with the _E engine (ckpt_segment_E /
//     replay_trajectory_cost). The two engines differ -> J(phi_ref) != 0 even at the
//     true phi. Here we generate the references with build_checkpoint_records (the EXACT
//     _E forward). With gen_ref_mode on, that forward dumps the _E state via
//     storeReferenceMap() at each obs trigger. Assimilation later replays the same _E
//     trajectory (same phi, same dt sequence, same obs grid) -> exact match -> J = 0.
// ===========================================================================
bool swSpE::generate_references_E(double* phi_map, size_t checkpoint_stride)
{
    const double dt_ref = dt_map_ad;
    if (dt_ref <= small_real) {
        logError() << "generate_references_E : dt_map_ad invalide (" << dt_ref << ")";
        return false;
    }

    // (Re)create reference/ + enable gen_ref_mode + reset n_ref_written.
    if (!begin_reference_generation()) return false;

    // Synthesize the observation-time grid: multiples of dt_ref inside the window
    // [t_state_assimilation, t_max] — identical to the cadence load_all_observations()
    // will re-read from the filenames, so assimilation lands on the SAME instants.
    obs_times_vec.clear();
    long k0 = (long)std::ceil((t_state_assimilation - small_real) / dt_ref);
    if (k0 < 0) k0 = 0;
    for (long k = k0; ; ++k) {
        const double tt = (double)k * dt_ref;
        if (tt > t_max + small_real) break;
        if (tt >= t_state_assimilation - small_real) obs_times_vec.push_back(tt);
    }
    if (obs_times_vec.empty()) {
        logError() << "generate_references_E : aucune date d'observation dans la fenetre ["
                   << t_state_assimilation << "s, " << t_max << "s] (dt_map_ad=" << dt_ref << ").";
        gen_ref_mode = false;
        return false;
    }
    current_obs_times  = obs_times_vec.data();
    current_obs_frames = nullptr;          // not read by build_checkpoint_records
    current_n_obs      = obs_times_vec.size();

    logInfo() << "generate_references_E : " << current_n_obs
              << " date(s) d'observation sur [" << t_state_assimilation << "s, " << t_max
              << "s] tous les " << dt_ref << "s ; forward _E (stride=" << checkpoint_stride << ").";

    // Freeze the normalize bounds for this phi (same as assimilation does before its tape),
    // so phi_correct(phi_map) is identical between generation and assimilation.
    update_phi_norm_bounds(phi_map);

    // Run the EXACT _E forward. gen_ref_mode -> storeReferenceMap() fires at each trigger.
    build_checkpoint_records(phi_map, checkpoint_stride);

    // Disable gen_ref_mode + sanity-check the number of references written.
    return finalize_reference_generation();
}


// ===========================================================================
// 1ter. run_forward_E
//   Lance UNIQUEMENT le forward _E (build_checkpoint_records, la meme boucle que celle
//   differenciee par Enzyme) sur [t_0, t_max] et ecrit les cartes ResultHydro dans
//   output/ a la cadence dt_map. Aucune assimilation, aucun adjoint.
// ===========================================================================
bool swSpE::run_forward_E(double* phi_map, size_t checkpoint_stride)
{
    const double dt_out = (dt_map > small_real) ? dt_map : dt_max;
    if (dt_out <= small_real) {
        logError() << "run_forward_E : dt_map invalide (" << dt_out << ")";
        return false;
    }

    // (Re)cree output/ et active le dump forward vers output/.
    const std::string out_dir = m_io->absolutePath("output");
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    state_dump_dir   = "output";
    forward_dump_mode = true;
    n_ref_written     = 0;

    // Grille de dump : multiples de dt_map sur [t_0, t_max].
    obs_times_vec.clear();
    long k0 = (long)std::ceil((t_0 - small_real) / dt_out);
    if (k0 < 0) k0 = 0;
    for (long k = k0; ; ++k) {
        const double tt = (double)k * dt_out;
        if (tt > t_max + small_real) break;
        if (tt >= t_0 - small_real) obs_times_vec.push_back(tt);
    }
    if (obs_times_vec.empty()) {
        logError() << "run_forward_E : aucune date de sortie sur [" << t_0 << "s, " << t_max
                   << "s] (dt_map=" << dt_out << ").";
        forward_dump_mode = false; state_dump_dir = "reference";
        return false;
    }
    current_obs_times  = obs_times_vec.data();
    current_obs_frames = nullptr;            // non lu par build_checkpoint_records
    current_n_obs      = obs_times_vec.size();

    logInfo() << "run_forward_E : forward _E sur [" << t_0 << "s, " << t_max << "s], "
              << current_n_obs << " cartes tous les " << dt_out << "s -> " << out_dir;

    // Fige les bornes normalize pour ce phi (coherent avec gen-ref / assimilation).
    update_phi_norm_bounds(phi_map);

    // LA boucle temporelle differenciee par Enzyme.
    build_checkpoint_records(phi_map, checkpoint_stride);

    forward_dump_mode = false;
    state_dump_dir    = "reference";
    logInfo() << "run_forward_E : " << n_ref_written << " carte(s) ResultHydro ecrite(s) dans "
              << out_dir << ".";
    return true;
}


// ===========================================================================
// 2.  ckpt_segment_E  (static — no hidden 'this')
//     This is the function Enzyme differentiates for each segment.
//
//     enzyme_dup  : solver, phi_map, h_in, qx_in, qy_in
//     enzyme_const: dt_seq, bc_flat, obs_trig, n_steps, edge_bnd_nb,
//                   obs_times, obs_frames, n_obs, t_start
// ===========================================================================
double swSpE::ckpt_segment_E(
    swSpE*        solver,
    double*      phi_map,
    const double* h_in,
    const double* qx_in,
    const double* qy_in,
    const double* dt_seq,
    const double* bc_flat,
    const int*    obs_trig,
    size_t        n_steps,
    size_t        ebnd_nb,
    const double* obs_times,
    const double* const* obs_frames,
    size_t        n_obs,
    double        t_start)
{
    // --- Load checkpoint state into solver ---
    // (Enzyme differentiates through these assignments: d_h_in accumulates here)
    for (size_t i = 0; i < solver->cell_nb; i++) {
        solver->water_depth[i]      = h_in[i];
        solver->unit_discharge_x[i] = qx_in[i];
        solver->unit_discharge_y[i] = qy_in[i];
        // CHANGEMENT DE VARIABLE DANS LA TAPE : phi_map est la variable d'optimisation
        // m (non bornee). phi_correct applique le mode choisi (none/clamp/tanh/normalize).
        // Enzyme derive dphi/dm automatiquement (pour normalize, lo/hi sont figes en amont
        // -> branche affine -> derivee 1/(hi-lo)), et computeEdgeFactors_E (juste apres)
        // voit le phi transforme.
        solver->phi[i]              = solver->phi_correct(phi_map[i]);
    }
    // FIX Bug 1 : ne PAS aliaser phi_ptr sur phi_map.
    // phi[] vient d'être chargé depuis phi_map dans la boucle ci-dessus.
    // Si on écrivait "solver->phi_ptr = phi_map", Enzyme compterait le gradient
    // de φ deux fois : une fois via l'argument explicite enzyme_dup phi_map/dJ_dPhi,
    // une fois via solver->phi_ptr (partie de enzyme_dup solver/adj).
    // En pointant sur solver->phi.data() on n'a qu'un seul chemin → pas de doublons.
    solver->phi_ptr = solver->phi.data();
    solver->t       = t_start;

    // FIX C : recalcul des facteurs 1/phi DANS la zone differenciee, a partir du
    // phi qu'on vient de charger (ligne ci-dessus). Avant ce fix, les facteurs
    // etaient figes (calcules a l'init / dans build_checkpoint_records, hors tape)
    // -> Enzyme les voyait comme des constantes et ratait d(1/phi)/dphi = -1/phi^2.
    // En les recalculant ici, phi -> coeff -> bilan -> cout est une chaine de
    // donnees taped -> le gradient capte enfin le levier dominant de la porosite.
    solver->computeEdgeFactors_E();

    // --- Enzyme activity analysis note ---
    // Enzyme performs interprocedural activity analysis: it only allocates tape
    // entries for values whose adjoints can be non-zero, i.e. values whose
    // computational path leads back to an active input (phi_map here).
    //
    // The following data in 'solver' is GEOMETRICALLY STATIC and never written
    // during the simulation:
    //   cell_z_b, cell_area, cell_centre_x/y, edge_int_norm_x/y,
    //   edge_bnd_norm_x/y, boussinesq, friction_cos/sin_alpha,
    //   friction_coeff_direction_1/2, cell_neigh_cell_idx
    //   (NB FIX C : edge_*_factor_* ne sont PLUS statiques — ils sont recalcules
    //    depuis phi par computeEdgeFactors_E() ci-dessus, donc ACTIFS pour Enzyme.)
    //
    // Because none of these derive from phi_map, Enzyme's forward activity
    // analysis marks them inactive → no shadow allocation, no tape entries.
    //
    // The data that IS active (derives from phi_map):
    //   water_depth, unit_discharge_x/y, phi[],
    //   and the intermediate Riemann states / fluxes written each step.
    //
    // The data we explicitly prevent from being tracked via enzyme_const args:
    //   dt_seq, bc_flat, obs_trig — prerecorded from the dry run.

    double cost = 0.0;

    for (size_t step = 0; step < n_steps; step++)
    {
        // ----- Prerecorded dt (enzyme_const) -----
        double cur_dt = dt_seq[step];
        solver->dt    = cur_dt;  // needed by MUSCL EVR reconstruction

        // ----- Prerecorded BC values (enzyme_const) -----
        // bc_val_interpolated is read by hypComputeFluxes_E for boundary edges.
        // Writing from an enzyme_const source: Enzyme won't differentiate through BC data.
        const double* bc_step = bc_flat + step * ebnd_nb;
        for (size_t j = 0; j < ebnd_nb; j++) {
            solver->bc_val_interpolated[j] = bc_step[j];
        }

        // ----- Physics (Enzyme WILL differentiate through all of this) -----
        solver->hypComputeLRStates_E();
        solver->hypComputeFluxes_E();
        solver->hypComputeBalance_E(cur_dt);

        // Dry-cell floor (Enzyme-safe — just sets to zero via conditional)
        for (size_t n = 0; n < solver->cell_nb; ++n) {
            if (solver->water_depth[n] < 0.) {
                solver->water_depth[n]      = 0.;
                solver->unit_discharge_x[n] = 0.;
                solver->unit_discharge_y[n] = 0.;
            }
        }

        // FIX cohérence checkpoint/replay : n'appliquer la friction QUE si
        // model_friction != "none", EXACTEMENT comme build_checkpoint_records (le dry-run
        // qui enregistre h_ckpt). Sinon sourceFrictionStrickler_E injecte une friction
        // fantôme (K=20 par défaut, coeffs non initialisés) que le dry-run, lui, ne met pas
        // -> les checkpoints h_ckpt[k] ne sont PAS sur la trajectoire de ckpt_segment_E ->
        // l'adjoint (qui restaure h_ckpt[k] à chaque segment) se linéarise au mauvais état ->
        // biais de gradient croissant avec le nombre de segments (Taylor plateau ∝ #segments).
        // model_friction est enzyme_const (ne dérive pas de phi) -> branche constante, tape OK.
        if (solver->model_friction != "none")
            solver->sourceFrictionStrickler_E(cur_dt);
        // Head loss singulier APRES la friction, IDENTIQUE au dry-run build_checkpoint_records.
        // model_singular_head_loss est enzyme_const -> branche constante, tape OK. Terme de
        // dissipation q*=exp(-AA*dt) : necessaire a la coherence physique ET stabilise
        // l'adjoint (sans lui, cas Toce : explosion adjointe -> gradient faux).
        if (solver->model_singular_head_loss != "none")
            solver->sourceSingularHeadLoss_E(cur_dt);
        solver->t += cur_dt;

        // ----- Observation cost -----
        // obs_type est enzyme_const (ne derive pas de phi) -> branche constante, tape OK
        // (meme mecanisme que le garde model_friction ci-dessus). Quand obs_type=="cells"
        // (defaut), le chemin SP/SP est strictement inchange.
        int oi = obs_trig[step];
        if (oi >= 0 && (size_t)oi < n_obs) {
            if (solver->obs_type == "gauges")
                cost += solver->compute_cost_gauges_ptr(obs_frames[oi]);
            else
                cost += solver->compute_cost_frame_ptr(obs_frames[oi]);
        }
    }

    return cost;
}


// ===========================================================================
// 3.  get_spatial_gradient_checkpointed
//     Outer driver: dry run → backward loop over segments.
// ===========================================================================
void swSpE::get_spatial_gradient_checkpointed(double* phi_map,
                                              double* dJ_dPhi,
                                              size_t  checkpoint_stride)
{
    const size_t N = cell_nb;

    // -----------------------------------------------------------------------
    // (a) Dry forward run — records ALL segments
    // -----------------------------------------------------------------------
    adjoint_truncation_seg_ = -1;  // reset from any previous run

    // Mode normalize : on (re)calcule lo/hi = min/max de m sur les cellules libres a
    // partir de la carte courante, UNE fois par iteree (= ici, ou grad() est appele),
    // puis on les fige. Toute la suite (build_checkpoint_records ci-dessous, la tape
    // ckpt_segment_E, et les replay_trajectory_cost de la recherche lineaire) utilise
    // ces bornes figees -> f et grad differencient la meme fonction (test de Taylor OK)
    // et Enzyme voit (m-lo)/(hi-lo) comme affine en m.
    update_phi_norm_bounds(phi_map);

    build_checkpoint_records(phi_map, checkpoint_stride);

    if (segs_record_.empty()) {
        std::cout << "[Checkpoint] WARNING: no segments recorded.\n";
        return;
    }

    // -----------------------------------------------------------------------
    // (b) Allocate adjoint solver (shadow copy of primal)
    // -----------------------------------------------------------------------
    swSpE* adj = new swSpE();
    adj->allocate_as_adjoint(this);
    adj->phi_ptr = dJ_dPhi;  // gradient accumulates directly into output array

    // Wire up observation pointers into the adjoint solver (enzyme_const data)
    adj->current_obs_times  = this->current_obs_times;
    adj->current_obs_frames = this->current_obs_frames;
    adj->current_n_obs      = this->current_n_obs;

    // -----------------------------------------------------------------------
    // (c) Adjoint state carried backward between segments:
    //     d_h[i] = ∂(J_{k+1}+...+J_K) / ∂water_depth_at_end_of_seg_k
    // -----------------------------------------------------------------------
    std::vector<double> d_h (N, 0.0);  // starts at zero (no "future" after last seg)
    std::vector<double> d_qx(N, 0.0);
    std::vector<double> d_qy(N, 0.0);

    // Adjoint of the INPUT checkpoint state (output of __enzyme_autodiff)
    std::vector<double> d_h_in (N, 0.0);
    std::vector<double> d_qx_in(N, 0.0);
    std::vector<double> d_qy_in(N, 0.0);

    // Zero the output gradient (will be accumulated across all segments)
    std::fill(dJ_dPhi, dJ_dPhi + N, 0.0);

    // -----------------------------------------------------------------------
    // (d) Backward loop: k = K-1 … first_seg_in_window
    // -----------------------------------------------------------------------
    const size_t K = segs_record_.size();

    // -----------------------------------------------------------------------
    // FENÊTRE D'ASSIMILATION DYNAMIQUE (windowed_assimilation)
    // -----------------------------------------------------------------------
    // Quand windowed_assimilation == true, l'adjoint ne remonte plus jusqu'au
    // segment 0 (t = t_0) mais s'arrête à first_seg_in_window, ce qui réalise une
    // rétropropagation sur la seule fenêtre [t_start_window, t_max].
    //
    // La taille de la fenêtre est pilotée par input.solv, de deux manières
    // mutuellement exclusives :
    //   - assimilation_window_duration : durée physique en secondes
    //                                    (t_start_window = t_max - duration) ;
    //   - assimilation_window_steps    : nombre de pas de temps comptés à rebours
    //                                    depuis t_max.
    // À défaut, on retombe sur t_state_assimilation.
    //
    // La troncature réutilise EXACTEMENT le mécanisme adjoint_truncation_seg_ :
    // en le positionnant à (first_seg_in_window - 1), replay_trajectory_cost
    // (utilisé pour le coût J ET le test de Taylor, donc par PSQN func()/grad())
    // redémarre au même segment → coût et gradient restent cohérents sur la
    // fenêtre, sans aucune autre modification.
    // -----------------------------------------------------------------------
    size_t first_seg_in_window = 0;  // défaut : remontée complète jusqu'au segment 0
    if (windowed_assimilation && K > 0)
    {
        if (assimilation_window_steps > 0)
        {
            // Mode "nombre de pas" : on cumule les pas à rebours jusqu'à atteindre
            // assimilation_window_steps. La troncature se fait à la frontière de
            // segment (granularité = checkpoint_stride pas) couvrant >= N pas.
            long acc = 0;
            first_seg_in_window = 0;
            for (long k = (long)K - 1; k >= 0; --k) {
                acc += (long)segs_record_[(size_t)k].n_steps;
                if (acc >= assimilation_window_steps) { first_seg_in_window = (size_t)k; break; }
            }
            std::cout << "[Window] Mode 'pas de temps' : fenetre = "
                      << assimilation_window_steps << " pas demandes -> "
                      << acc << " pas effectifs (frontiere de segment).\n";
        }
        else if (assimilation_window_duration > 0.)
        {
            // Mode "durée" : t_start_window = t_max - duration. On choisit le segment
            // qui CONTIENT t_start_window (le dernier dont t_start <= t_start_window),
            // afin que la fenetre couvre bien tout [t_start_window, t_max].
            const double t_start_window = t_max - assimilation_window_duration;
            first_seg_in_window = 0;
            for (size_t k = 0; k < K; ++k) {
                if (segs_record_[k].t_start <= t_start_window + 1e-9) first_seg_in_window = k;
                else break;
            }
            if (!minim_quiet_logs)
            std::cout << "[Window] Mode 'duree' : fenetre = "
                      << assimilation_window_duration << "s -> t_start_window = "
                      << t_start_window << "s.\n";
        }
        else if (t_state_assimilation > t_0)
        {
            // Repli : on utilise t_state_assimilation comme borne de fenetre.
            const double t_start_window = t_state_assimilation;
            first_seg_in_window = 0;
            for (size_t k = 0; k < K; ++k) {
                if (segs_record_[k].t_start <= t_start_window + 1e-9) first_seg_in_window = k;
                else break;
            }
            if (!minim_quiet_logs)
            std::cout << "[Window] Repli sur t_state_assimilation = "
                      << t_start_window << "s.\n";
        }

        if (first_seg_in_window > 0) {
            // Aligne la machinerie de troncature : replay/Taylor redemarrent ici.
            adjoint_truncation_seg_ = (int)first_seg_in_window - 1;
            if (!minim_quiet_logs) {
            std::cout << "[Window] *** FENETRE D'ASSIMILATION ACTIVE ***\n"
                      << "         Adjoint restreint aux segments " << first_seg_in_window
                      << " ... " << (K - 1) << " (sur " << K << ").\n"
                      << "         Retropropagation : t_max=" << t_max
                      << "s -> t_start_window=" << segs_record_[first_seg_in_window].t_start
                      << "s (au lieu de t_0=" << t_0 << "s).\n";
            logInfo() << "[Window] Fenetre d'assimilation active : segments "
                      << first_seg_in_window << "..." << (K - 1)
                      << ", t in [" << segs_record_[first_seg_in_window].t_start
                      << "s, " << t_max << "s].";
            }
        } else {
            std::cout << "[Window] windowed_assimilation actif mais fenetre >= trajectoire "
                      << "complete : remontee jusqu'au segment 0 (comportement standard).\n";
        }
    }

    // -----------------------------------------------------------------------
    // FIX L-BFGS-B (gel) : si la fenetre de troncature est FIGEE (forced_first_seg_>=0,
    // pose par le solveur L-BFGS-B au segment d'explosion du point de depart), on force
    // la borne de fenetre a cette valeur pour TOUTES les evaluations. La troncature
    // dynamique d'explosion (plus bas) est desactivee dans ce mode -> J calcule sur la
    // MEME fenetre a chaque appel -> line-search coherente, plus de gel.
    // -----------------------------------------------------------------------
    if (forced_first_seg_ >= 0) {
        first_seg_in_window = (size_t)std::min<long>((long)forced_first_seg_, (long)K - 1);
        adjoint_truncation_seg_ = (int)first_seg_in_window - 1;
    }

    // -----------------------------------------------------------------------
    // FIX Bug 4 (revised) — "watermark" save for explosion handling.
    //
    // Problem: the adjoint explosion is GRADUAL.  Between the last milestone
    // (e.g. k=150, ||d_h||=191) and the detected explosion (k=112, ||d_h||=1e11),
    // there are ~37 unlisted segments where d_h grows exponentially but stays
    // below the detection threshold (1e8 × ref).  Their __enzyme_autodiff calls
    // accumulate LARGE, CORRUPT contributions into dJ_dPhi even though each
    // individual segment doesn't trigger the check.
    //
    // The dJ_save saved *just before* the exploded segment k=112 already
    // includes these corrupt contributions from segs 149→113.
    //
    // Solution: maintain a "watermark" dJ_dPhi that is refreshed at every
    // milestone (k % 50 == 0, obs segment, or k == K-1).  When explosion is
    // detected, restore from this watermark (which was last updated at k=150
    // with clean adjoint state ||d_h||=191), not from the per-call save.
    //
    // The watermark is updated AFTER the backward pass of each milestone
    // segment, so it includes that segment's clean contribution.  When restored,
    // adjoint_truncation_seg_ is set to (watermark_seg - 1) so that
    // replay_trajectory_cost starts from segs_record_[watermark_seg].h_ckpt,
    // exactly mirroring what Enzyme accumulated.
    // -----------------------------------------------------------------------
    std::vector<double> dJ_watermark(N, 0.0);  // latest "safe" dJ_dPhi snapshot
    int watermark_seg = (int)K - 1;            // segment index of last save (non-tiled)
    int watermark_kept_start = (int)K;         // 1er segment conserve si un pavé explose (tiled)

    // Static across gradient calls — persists between segments within one run.
    static double norm_dh_at_last_obs = -1.0;
    norm_dh_at_last_obs = -1.0;  // reset for this gradient computation
    double dh_peak = 0.0;        // pic de ||d_h|| sur le backward (diagnostic descente)

    // -----------------------------------------------------------------------
    // PAVAGE DE LA FENETRE (n_assimilation_windows > 1)
    // -----------------------------------------------------------------------
    // On decoupe la fenetre [first_seg_in_window, K-1] en n_assimilation_windows
    // pavés contigus de tailles (en nombre de SEGMENTS) aussi egales que possible
    // — choix le plus stable, car la longueur d'un adjoint ∝ nombre de segments,
    // donc tous les pavés partagent le meme horizon de stabilite.
    //
    // A chaque frontiere de pavé (en remontant) on remet l'etat adjoint d_h=0
    // (lambda=0 a la fin du pavé) : chaque pavé devient un adjoint court et borne,
    // et dJ_dPhi accumule la SOMME des gradients partiels (entree de pavé figee).
    //
    // assim_window_tiles_ (membre) memorise les segments de DEBUT de chaque pavé ;
    // replay_trajectory_cost l'utilise pour calculer le cout comme somme des couts
    // par pavé depuis checkpoints figes -> coherence func()/grad()/Taylor.
    // -----------------------------------------------------------------------
    assim_window_tiles_.clear();
    std::vector<char> tile_top(K, 0);   // tile_top[k]=1 => k est le dernier seg d'un pavé
    const bool tiled = (n_assimilation_windows > 1) && (first_seg_in_window < K);
    if (tiled) {
        const size_t a = first_seg_in_window;
        const size_t M = K - a;                       // nb de segments dans la fenetre
        size_t n = (size_t)n_assimilation_windows;
        if (n > M) n = M;                             // au plus 1 pavé par segment
        for (size_t tdx = 0; tdx < n; ++tdx) {
            const size_t s_t = a + (tdx * M) / n;         // debut du pavé tdx (division entiere)
            const size_t e_t = a + ((tdx + 1) * M) / n - 1; // dernier seg du pavé tdx
            assim_window_tiles_.push_back(s_t);
            tile_top[e_t] = 1;                            // reset d_h avant ce seg (en remontant)
        }
        if (!minim_quiet_logs) {
        std::cout << "[Window] PAVAGE : " << n << " pavés sur la fenetre [seg "
                  << a << ", seg " << (K - 1) << "] (~" << (M / n)
                  << " segments/pavé). Adjoint borne pavé par pavé.\n";
        logInfo() << "[Window] Pavage actif : " << n << " fenetres d'assimilation sur"
                  << " les segments " << a << ".." << (K - 1) << ".";
        }
    }

    for (int k = (int)K - 1; k >= (int)first_seg_in_window; k--)
    {
        const SegPhysData& seg = segs_record_[(size_t)k];

        // ---- Frontiere de pavé : remettre l'adjoint a zero (lambda=0 en fin de pavé) ----
        // Chaque pavé est un adjoint court et borne ; dJ_dPhi accumule la somme des
        // gradients partiels (etat d'entree du pavé traite comme fige). On prend ici un
        // snapshot "propre" : les pavés plus recents (deja traites) sont sains -> sert de
        // point de restauration si le pavé courant venait a exploser.
        if (tiled && tile_top[(size_t)k]) {
            std::fill(d_h.begin(),  d_h.end(),  0.0);
            std::fill(d_qx.begin(), d_qx.end(), 0.0);
            std::fill(d_qy.begin(), d_qy.end(), 0.0);
            std::copy(dJ_dPhi, dJ_dPhi + N, dJ_watermark.begin());
            watermark_kept_start = k + 1;   // debut du 1er pavé conserve si ce pavé explose
        }

        // ---- Seed adjoint solver's state with d_h/d_qx/d_qy from seg k+1 ----
        // These represent ∂(J_{k+1}+...)/∂(final state of seg k).
        // Enzyme will propagate them backward through seg k's steps.
        for (size_t i = 0; i < N; i++) {
            adj->water_depth[i]      = d_h[i];
            adj->unit_discharge_x[i] = d_qx[i];
            adj->unit_discharge_y[i] = d_qy[i];
        }
        // phi adjoint in adj->phi is already zero from allocate_as_adjoint;
        // Enzyme will ADD the gradient of seg k into adj->phi (= dJ_dPhi via phi_ptr).

        // ---- Reset d_h_in to zero: they'll accumulate adjoint of h_start of seg k ----
        std::fill(d_h_in.begin(),  d_h_in.end(),  0.0);
        std::fill(d_qx_in.begin(), d_qx_in.end(), 0.0);
        std::fill(d_qy_in.begin(), d_qy_in.end(), 0.0);

        // ---- Restore primal solver to checkpoint state ----
        // (Enzyme will replay from here during its internal forward pass)
        for (size_t i = 0; i < N; i++) {
            water_depth[i]      = seg.h_ckpt[i];
            unit_discharge_x[i] = seg.qx_ckpt[i];
            unit_discharge_y[i] = seg.qy_ckpt[i];
            phi[i]              = phi_map[i];
        }
        // FIX Bug 1 (primal side) : même raison que dans ckpt_segment_E.
        // phi[] vient d'être rechargé depuis phi_map; pointer sur phi.data()
        // évite l'aliasing avec l'argument enzyme_dup phi_map/dJ_dPhi.
        phi_ptr = phi.data();
        this->t = seg.t_start;

        // ---- Differentiate segment k ----
        // enzyme_dup  : primal/adjoint solver, phi_map/dJ_dPhi,
        //               h_ckpt/d_h_in, qx_ckpt/d_qx_in, qy_ckpt/d_qy_in
        // enzyme_const: all prerecorded data + topology scalars

        // FIX Bug 2 : sauvegarder dJ_dPhi AVANT l'appel Enzyme (rollback du segment k).
        std::vector<double> dJ_save(dJ_dPhi, dJ_dPhi + N);

        __enzyme_autodiff(
            (void*)swSpE::ckpt_segment_E,
            // solver (contains working state)
            enzyme_dup,   this,  adj,
            // phi (the parameter we optimise)
            enzyme_dup,   phi_map, dJ_dPhi,
            // input state at start of segment
            enzyme_dup,   seg.h_ckpt.data(),  d_h_in.data(),
            enzyme_dup,   seg.qx_ckpt.data(), d_qx_in.data(),
            enzyme_dup,   seg.qy_ckpt.data(), d_qy_in.data(),
            // prerecorded (enzyme_const) — Enzyme never tracks adjoints through these
            enzyme_const, seg.dt_seq.data(),
            enzyme_const, seg.bc_flat.data(),
            enzyme_const, seg.obs_trig.data(),
            enzyme_const, seg.n_steps,
            enzyme_const, this->edge_bnd_nb,
            enzyme_const, this->current_obs_times,
            enzyme_const, this->current_obs_frames,
            enzyme_const, this->current_n_obs,
            enzyme_const, seg.t_start
        );

        // ---- Harvest adjoint of input state → seed for segment k-1 ----
        // d_h_in[i] = ∂(J_k + J_{k+1} + ...)/∂h_start_k[i]
        //           = adjoint of the OUTPUT state of segment k-1
        for (size_t i = 0; i < N; i++) {
            d_h[i]  = d_h_in[i];
            d_qx[i] = d_qx_in[i];
            d_qy[i] = d_qy_in[i];
        }

        // Reset adjoint phi in adj (it was accumulated into dJ_dPhi via phi_ptr).
        // Must be zero for the next segment call so gradients don't double-count.
        std::fill(adj->phi.begin(), adj->phi.end(), 0.0);

        // --- Norme L2 de l'état adjoint ---
        double norm_dh = 0.0;
        for (size_t i = 0; i < N; i++) norm_dh += d_h[i] * d_h[i];
        norm_dh = std::sqrt(norm_dh);
        if (norm_dh > dh_peak) dh_peak = norm_dh;

        // --- Mise à jour de la référence (premier segment non-nul) ---
        if (norm_dh_at_last_obs < 0.0 && norm_dh > 0.0) {
            norm_dh_at_last_obs = norm_dh;
        }

        // --- Comptage observations dans ce segment ---
        int n_obs_in_seg = 0;
        for (size_t s = 0; s < seg.n_steps; s++)
            if (seg.obs_trig[s] >= 0) n_obs_in_seg++;

        bool has_obs     = (n_obs_in_seg > 0);
        bool is_milestone = (k % 50 == 0) || (k == (int)K-1) || (k == 0);

        // --- Watermark + détection d'explosion (arret precoce de l'adjoint) ---
        bool adjoint_blown_up = false;

        if (!tiled) {
        // ===== FENETRE UNIQUE / pas de pavage : comportement historique =====
        if (is_milestone || has_obs) {
            std::copy(dJ_dPhi, dJ_dPhi + N, dJ_watermark.begin());
            watermark_seg = k;
        }
        if (forced_first_seg_ < 0 && norm_dh_at_last_obs > 0.0 && norm_dh > 1e8 * norm_dh_at_last_obs) {
            adjoint_blown_up = true;
            // On saute une ligne pour ne pas écraser la ligne d'avancement dynamique
            std::cout << "\n[Checkpoint] *** EXPLOSION ADJOINT au segment " << k
                    << " : ||d_h||=" << norm_dh
                    << " >> " << norm_dh_at_last_obs << " (ref)\n"
                    << "             Dernier jalon propre : segment " << watermark_seg << ".\n"
                    << "             Gradients restreints aux segments "
                    << watermark_seg << " … " << (K-1) << ".\n"
                    << "             Cause : instabilité de l'adjoint de St-Venant sur\n"
                    << "             longue fenêtre (ondes amplifiées en sens inverse du temps).\n"
                    << "             Remède : activer le PAVAGE (windowed_assimilation=1 +\n"
                    << "             n_assimilation_windows>1) ou reduire la fenetre.\n";
            logWarn() << "[Checkpoint] *** ARRET PRECOCE DE L'ADJOINT : explosion au segment "
                    << k << " (||d_h||=" << norm_dh << " >> ref=" << norm_dh_at_last_obs
                    << "). Gradient restreint aux segments " << watermark_seg << ".." << (K-1)
                    << ". Activer le pavage (n_assimilation_windows) ou reduire la fenetre.";
            std::copy(dJ_watermark.begin(), dJ_watermark.end(), dJ_dPhi);
            adjoint_truncation_seg_ = watermark_seg - 1;
            break;
        }
    } else {
        // ===== MODE PAVE : une explosion signifie un pavé encore trop long =====
        if (forced_first_seg_ < 0 && norm_dh_at_last_obs > 0.0 && norm_dh > 1e8 * norm_dh_at_last_obs) {
            adjoint_blown_up = true;
            // On saute une ligne pour ne pas écraser la ligne d'avancement dynamique
            std::cout << "\n[Checkpoint] *** EXPLOSION ADJOINT dans un pavé (segment " << k
                    << ") : ||d_h||=" << norm_dh
                    << " >> " << norm_dh_at_last_obs << " (ref)\n"
                    << "             Pavés conserves : segments >= " << watermark_kept_start
                    << " (les plus proches de t_max).\n"
                    << "             Remède : AUGMENTER n_assimilation_windows"
                    << " (pavés plus courts => adjoint plus stable).\n";
            logWarn() << "[Checkpoint] *** ARRET PRECOCE (pavé instable) au segment " << k
                    << " (||d_h||=" << norm_dh << " >> ref=" << norm_dh_at_last_obs
                    << "). Pavés conserves >= seg " << watermark_kept_start
                    << ". Augmenter n_assimilation_windows.";
            std::copy(dJ_watermark.begin(), dJ_watermark.end(), dJ_dPhi);
            adjoint_truncation_seg_ = watermark_kept_start - 1;
            while (!assim_window_tiles_.empty() &&
                assim_window_tiles_.front() < (size_t)watermark_kept_start)
                assim_window_tiles_.erase(assim_window_tiles_.begin());
            break;
        }
    }

    // ===== AFFICHAGE DYNAMIQUE SUR UNE SEULE LIGNE =====
    // Note : J'ai enlevé la condition restrictive pour que tu voies l'avancement défiler fluidement,
    // mais tu peux la remettre si tu veux que ça se rafraîchisse moins souvent.
    if (has_obs || is_milestone || (tiled && tile_top[(size_t)k])) {
        std::cout << "\r[Checkpoint] Segment " << k << "/" << (K-1)
                << "  n_obs=" << n_obs_in_seg << "   " << std::flush;
    }

        // seulement affichage :
        /*
        if (!tiled) {
            // ===== FENETRE UNIQUE / pas de pavage : comportement historique =====
            // Watermark rafraichi aux jalons (milestone + obs) ; a ce stade dJ_dPhi
            // contient les contributions des segments K-1 → k (stables jusqu'ici).
            if (is_milestone || has_obs) {
                std::copy(dJ_dPhi, dJ_dPhi + N, dJ_watermark.begin());
                watermark_seg = k;
            }
            if (forced_first_seg_ < 0 && norm_dh_at_last_obs > 0.0 && norm_dh > 1e8 * norm_dh_at_last_obs) {
                adjoint_blown_up = true;
                std::cout << "[Checkpoint] *** EXPLOSION ADJOINT au segment " << k
                          << " : ||d_h||=" << norm_dh
                          << " >> " << norm_dh_at_last_obs << " (ref)\n"
                          << "             Dernier jalon propre : segment " << watermark_seg << ".\n"
                          << "             Gradients restreints aux segments "
                          << watermark_seg << " … " << (K-1) << ".\n"
                          << "             Cause : instabilité de l'adjoint de St-Venant sur\n"
                          << "             longue fenêtre (ondes amplifiées en sens inverse du temps).\n"
                          << "             Remède : activer le PAVAGE (windowed_assimilation=1 +\n"
                          << "             n_assimilation_windows>1) ou reduire la fenetre.\n";
                logWarn() << "[Checkpoint] *** ARRET PRECOCE DE L'ADJOINT : explosion au segment "
                          << k << " (||d_h||=" << norm_dh << " >> ref=" << norm_dh_at_last_obs
                          << "). Gradient restreint aux segments " << watermark_seg << ".." << (K-1)
                          << ". Activer le pavage (n_assimilation_windows) ou reduire la fenetre.";
                std::copy(dJ_watermark.begin(), dJ_watermark.end(), dJ_dPhi);
                adjoint_truncation_seg_ = watermark_seg - 1;
                break;
            }
        } else {
            // ===== MODE PAVE : une explosion signifie un pavé encore trop long =====
            // Le watermark propre a deja ete pris a la frontiere du pavé courant (en
            // haut de boucle). On restaure les pavés deja propres (plus proches de
            // t_max), on coupe ici (arret precoce au grain du pavé) et on tronque la
            // liste des pavés pour garder la coherence du cout (replay_trajectory_cost).
            if (forced_first_seg_ < 0 && norm_dh_at_last_obs > 0.0 && norm_dh > 1e8 * norm_dh_at_last_obs) {
                adjoint_blown_up = true;
                std::cout << "[Checkpoint] *** EXPLOSION ADJOINT dans un pavé (segment " << k
                          << ") : ||d_h||=" << norm_dh
                          << " >> " << norm_dh_at_last_obs << " (ref)\n"
                          << "             Pavés conserves : segments >= " << watermark_kept_start
                          << " (les plus proches de t_max).\n"
                          << "             Remède : AUGMENTER n_assimilation_windows"
                          << " (pavés plus courts => adjoint plus stable).\n";
                logWarn() << "[Checkpoint] *** ARRET PRECOCE (pavé instable) au segment " << k
                          << " (||d_h||=" << norm_dh << " >> ref=" << norm_dh_at_last_obs
                          << "). Pavés conserves >= seg " << watermark_kept_start
                          << ". Augmenter n_assimilation_windows.";
                std::copy(dJ_watermark.begin(), dJ_watermark.end(), dJ_dPhi);
                adjoint_truncation_seg_ = watermark_kept_start - 1;
                while (!assim_window_tiles_.empty() &&
                       assim_window_tiles_.front() < (size_t)watermark_kept_start)
                    assim_window_tiles_.erase(assim_window_tiles_.begin());
                break;
            }
        }

        // Affichage partiel (milestone, segs avec obs, frontieres de pavé)
        if (has_obs || is_milestone || adjoint_blown_up || (tiled && tile_top[(size_t)k])) {
            std::cout << "[Checkpoint] Segment " << k << "/" << (K-1)
                      << "  n_steps=" << seg.n_steps
                      << "  ||d_h||=" << norm_dh
                      << "  n_obs=" << n_obs_in_seg
                      << (has_obs ? "  *** OBS ***" : "")
                      << "\n";
        }*/
    }

    
    {
        size_t n_masked = 0, n_wall = 0;
        for (size_t i = 0; i < N; i++) {
            // (1) Assainissement TOUJOURS actif : NaN/inf -> 0 (norme/Taylor propres).
            if (std::isnan(dJ_dPhi[i]) || std::isinf(dJ_dPhi[i])) {
                dJ_dPhi[i] = 0.0;
                n_masked++;
                continue;
            }
            // (2) MASQUE MURS (phi<=phi_min) : OPTIONNEL, pilote par input.solv
            // 'mask_wall_gradient' (defaut FALSE = desactive). Sert a eviter la singularite
            // 1/phi en phi=0 ; desactive, les cellules a phi tres bas gardent leur gradient
            // (utile pour laisser le barrage se former), au prix d'un risque numerique pres
            // de phi=0 (rattrape par l'assainissement NaN/inf ci-dessus).
            if (mask_wall_gradient) {
                double phi_t = phi_correct(phi_map[i]);
                if (!(phi_t > phi_min)) { dJ_dPhi[i] = 0.0; n_masked++; n_wall++; }
            }
        }
        if (n_masked)
            std::cout << "  [Masque gradient] " << n_masked << " cellules forcees a 0"
                      << " (dont " << n_wall << " murs phi<=phi_min ; mask_wall_gradient="
                      << (mask_wall_gradient ? "1" : "0") << ")\n";
    }

    // -----------------------------------------------------------------------
    // (e) Bilan final du gradient
    // -----------------------------------------------------------------------
    // Norme L2 et distribution du gradient
    double norm_grad = 0.0;
    size_t n_nonzero = 0;
    double grad_max = 0.0;
    size_t i_max = 0;
    for (size_t i = 0; i < N; i++) {
        double g = dJ_dPhi[i];
        norm_grad += g * g;
        if (g != 0.0) n_nonzero++;
        if (std::abs(g) > std::abs(grad_max)) { grad_max = g; i_max = i; }
    }
    norm_grad = std::sqrt(norm_grad);

    // Expose les diagnostics de cette evaluation (suivi de descente, ecrit dans le CSV).
    diag_grad_l2     = norm_grad;
    diag_grad_max    = grad_max;
    diag_grad_argmax = (long)i_max;
    diag_dh_peak     = dh_peak;

    if (!minim_quiet_logs) {
    std::cout << "=== GRADIENT CHECKPOINTED ===" << std::endl;
    std::cout << "  ||dJ/dPhi||_2 = " << norm_grad << "\n";
    std::cout << "  Cellules non-nulles : " << n_nonzero << " / " << N << "\n";
    std::cout << "  Valeur max : dJ/dphi[" << i_max << "] = " << grad_max << "\n";
    // Premières cellules non-nulles (pour vérification)
    size_t printed = 0;
    for (size_t i = 0; i < N && printed < 5; i++) {
        if (dJ_dPhi[i] != 0.0) {
            std::cout << "  dJ/dphi[" << i << "] = " << dJ_dPhi[i] << "\n";
            printed++;
        }
    }
    }
    if (n_nonzero == 0) {
        std::cout << "  ATTENTION : gradient entierement nul.\n"
                  << "  -> Verifier que les observations se declenchent (bilan forward ci-dessus).\n"
                  << "  -> Les obs à t=0 ne peuvent pas déclencher (dt_to_obs=0), c'est normal.\n";
    }
    if (norm_grad > 1e10) {
        std::cout << "  ATTENTION : norme du gradient tres grande (" << norm_grad << ").\n"
                  << "  L'adjoint a peut-etre explosé (instabilité sur longue fenêtre).\n"
                  << "  -> Réduire la fenêtre d'assimilation (ex: ne prendre que les\n"
                  << "     observations dans la derniere heure avant t_max).\n";
    }

    // Lancer le test de Taylor si le gradient est non-nul et de magnitude raisonnable
    // (désactivé pendant la minimisation via taylor_test_after_gradient = false).
    if (n_nonzero > 0 && taylor_test_after_gradient) {
        // n_top_cells=1 : test only the top-gradient cell + directional test.
        // For large meshes (Gardon: 583 segments) each replay_trajectory_cost
        // call takes ~180s → keep n_top_cells small for validation.
        // Increase to 3+ for detailed debugging on small test cases.
        run_taylor_test_checkpointed(phi_map, dJ_dPhi, /*n_top_cells=*/1);
    }

    delete adj;
}


// ===========================================================================
// 4.  replay_trajectory_cost
//     Évalue J en rejouant la TRAJECTOIRE FIXE enregistrée dans segs_record_,
//     de manière SÉQUENTIELLE (état transmis d'un segment au suivant).
//
//     Pourquoi séquentiel et non par checkpoint fixe ?
//     ─────────────────────────────────────────────────
//     Le backward checkpointé (get_spatial_gradient_checkpointed) calcule le
//     gradient TOTAL ∂J/∂φ via la règle de la chaîne sur TOUS les segments :
//
//       grad_total = grad_direct_seg_k
//                  + Σ_{k'<k} ∂J_k/∂h_start_k * ∂h_end_seg_k'/∂φ   (indirect)
//
//     Si replay_trajectory_cost utilisait des checkpoints FIXES (h_ckpt_k
//     indépendants de φ), perturber φ dans le segment k' ne changerait PAS
//     h_ckpt_k (fixé une fois pour toutes) → J_replay ne capturerait que la
//     contribution DIRECTE de φ dans chaque segment.  Le test de Taylor
//     donnerait alors ratio = grad_direct / grad_total < 1 (≈ 0.84 observé).
//
//     Avec la version SÉQUENTIELLE, l'état à la fin du segment k' est
//     retransmis comme état initial du segment k'+1 (comme dans la vraie
//     simulation).  Ainsi J_replay_seq capture exactement la même fonction
//     que le backward Enzyme → le test de Taylor donne ratio → 1.
//
//     Note : les dt/BC/obs_trig restent ceux du dry-run (fixes), donc J reste
//     lisse en φ (pas de sauts de CFL). Seuls h/qx/qy se propagent.
// ===========================================================================
double swSpE::replay_trajectory_cost(double* phi_map)
{
    if (segs_record_.empty()) return 0.0;

    // -----------------------------------------------------------------------
    // MODE PAVE (assim_window_tiles_ non vide) : cout = SOMME des couts par pavé,
    // chacun rejoue depuis son checkpoint d'entree FIGE (etat independant de la
    // perturbation de phi -> chemin inter-pavé coupe, exactement comme l'adjoint
    // pavé qui remet d_h=0 a chaque frontiere). Les frontieres sont celles
    // calculees par get_spatial_gradient_checkpointed -> func()/grad()/Taylor
    // differencient la MEME fonction. Avec 1 pavé, ce bloc redonne le replay
    // continu classique (mais on garde le chemin legacy ci-dessous pour ce cas).
    // -----------------------------------------------------------------------
    if (!assim_window_tiles_.empty()) {
        const size_t Ksegs = segs_record_.size();
        double total_cost = 0.0;
        for (size_t ti = 0; ti < assim_window_tiles_.size(); ++ti) {
            const size_t s_t = assim_window_tiles_[ti];
            const size_t e_t = (ti + 1 < assim_window_tiles_.size())
                                 ? assim_window_tiles_[ti + 1] - 1
                                 : Ksegs - 1;
            // Etat d'entree FIGE du pavé (checkpoint du dry-run).
            const SegPhysData& seg_start = segs_record_[s_t];
            for (size_t i = 0; i < cell_nb; i++) {
                water_depth[i]      = seg_start.h_ckpt[i];
                unit_discharge_x[i] = seg_start.qx_ckpt[i];
                unit_discharge_y[i] = seg_start.qy_ckpt[i];
                phi[i]              = phi_map[i];
            }
            phi_ptr = phi.data();
            this->t = seg_start.t_start;
            // Rejoue les segments du pavé en propageant l'etat A L'INTERIEUR du pavé
            // (water_depth.data() comme h_in => no-op de chargement, etat conserve).
            for (size_t si = s_t; si <= e_t; si++) {
                const auto& seg = segs_record_[si];
                total_cost += ckpt_segment_E(
                    this, phi_map,
                    water_depth.data(), unit_discharge_x.data(), unit_discharge_y.data(),
                    seg.dt_seq.data(), seg.bc_flat.data(), seg.obs_trig.data(),
                    seg.n_steps, this->edge_bnd_nb,
                    this->current_obs_times, this->current_obs_frames, this->current_n_obs,
                    seg.t_start);
            }
        }
        return total_cost;
    }

    // -----------------------------------------------------------------------
    // MODE LEGACY (replay continu) — fenetre unique / pas de pavage.
    // Déterminer le segment de départ.
    //
    // Si l'adjoint a explosé au segment k (adjoint_truncation_seg_ = k),
    // Enzyme n'a calculé le gradient que pour les segments k+1 … K-1, en
    // traitant l'état au checkpoint k+1 comme FIXE (d_h zeroed at segment k).
    //
    // Pour que le FD mesure la MÊME fonction de phi qu'Enzyme, on démarre
    // depuis ce checkpoint fixe et on rejoue uniquement les segments k+1 … K-1.
    // Ainsi phi[0…k] n'a aucun effet sur J_replay (le chemin indirect est coupé)
    // exactement comme dans le gradient Enzyme.
    //
    // Sans explosion (adjoint_truncation_seg_ = -1), on repart du segment 0
    // comme avant → comportement inchangé (cas E3T).
    // -----------------------------------------------------------------------
    size_t start_seg = 0;
    if (adjoint_truncation_seg_ >= 0 &&
        (size_t)adjoint_truncation_seg_ + 1 < segs_record_.size())
    {
        start_seg = (size_t)adjoint_truncation_seg_ + 1;
        // Message une seule fois par appel pour le diagnostic
        static bool first_truncated_replay = true;
        if (first_truncated_replay) {
            std::cout << "[Replay] Troncature active : démarrage au segment "
                      << start_seg << " (explosion adjoint au seg "
                      << adjoint_truncation_seg_ << ").\n"
                      << "         J_replay ne compte que les segments "
                      << start_seg << " … " << (segs_record_.size()-1)
                      << " — cohérent avec le gradient Enzyme tronqué.\n";
            first_truncated_replay = false;
        }
    }

    // --- Charger l'état initial depuis le checkpoint de départ (FIXE en phi) ---
    const SegPhysData& seg_start = segs_record_[start_seg];
    for (size_t i = 0; i < cell_nb; i++) {
        water_depth[i]      = seg_start.h_ckpt[i];
        unit_discharge_x[i] = seg_start.qx_ckpt[i];
        unit_discharge_y[i] = seg_start.qy_ckpt[i];
        phi[i]              = phi_map[i];
    }
    phi_ptr = phi.data();
    this->t = seg_start.t_start;

    double total_cost = 0.0;

    // --- Rejouer chaque segment en transmettant l'état courant ---
    // L'astuce : on passe water_depth.data() comme h_in.
    // ckpt_segment_E fait d'abord "solver->water_depth[i] = h_in[i]" : si
    // h_in pointe sur le même vecteur, c'est un no-op et l'état courant est
    // conservé.  Après les pas de physique, water_depth contient l'état de
    // FIN du segment → transmis automatiquement au segment suivant.
    for (size_t si = start_seg; si < segs_record_.size(); si++) {
        const auto& seg = segs_record_[si];
        total_cost += ckpt_segment_E(
            this,
            phi_map,
            water_depth.data(),       // état courant (transmis entre segments)
            unit_discharge_x.data(),  // idem
            unit_discharge_y.data(),  // idem
            seg.dt_seq.data(),
            seg.bc_flat.data(),
            seg.obs_trig.data(),
            seg.n_steps,
            this->edge_bnd_nb,
            this->current_obs_times,
            this->current_obs_frames,
            this->current_n_obs,
            seg.t_start
        );
        // À ce stade, water_depth/unit_discharge sont l'état FIN du segment.
        // phi est rechargé depuis phi_map au début de ckpt_segment_E → OK.
    }
    return total_cost;
}


// ===========================================================================
// 5.  run_taylor_test_checkpointed
//
//     POURQUOI utiliser replay_trajectory_cost et non execute_E ?
//     ─────────────────────────────────────────────────────────────
//     execute_E recalcule la CFL à chaque pas de temps. Quand on perturbe
//     phi d'une quantité ε, la CFL change à chaque pas → les dt divergent
//     entre les trajectoires phi et phi+ε·eᵢ → les observations déclenchent
//     à des instants différents → J(phi+ε) est discontinu en phi, même pour
//     ε minuscule → le test de Taylor donne des ratios chaotiques.
//
//     replay_trajectory_cost utilise les dt/BC/obs_trig préenregistrés
//     (fixes, indépendants de phi) → J est lisse → test de Taylor propre.
//
//     Théorie : si ∂J/∂φᵢ est exact,
//       r(ε) = [J_replay(φ+ε·eᵢ) - J_replay(φ)] / (ε·∂J/∂φᵢ)  →  1  quand ε→0
//     Convergence attendue d'ordre 1 : |r-1| ∝ ε (effets non-linéaires).
// ===========================================================================
void swSpE::run_taylor_test_checkpointed(double* phi_base,
                                         double* dJ_dPhi,
                                         size_t  n_top_cells)
{
    const size_t N = cell_nb;

    std::cout << "\n========================================\n"
              << " TEST DE TAYLOR — gradient checkpointé \n"
              << "========================================\n"
              << " r(ε) = [J_traj(φ+ε·eᵢ) - J_traj(φ)] / (ε·∂J/∂φᵢ)  → 1\n"
              << " J_traj = trajectoire FIXE (dt/obs préenregistrés) — cohérent avec Enzyme\n\n";

    // ------------------------------------------------------------------
    // 1. Sélection des n_top_cells cellules par |gradient| décroissant
    //    + n_top_cells cellules supplémentaires à magnitude moyenne
    //    (pour détecter si seules les cellules "max" sont correctes)
    // ------------------------------------------------------------------
    std::vector<std::pair<double, size_t>> sorted_grad;
    sorted_grad.reserve(N);
    for (size_t i = 0; i < N; i++) {
        if (dJ_dPhi[i] != 0.0)
            sorted_grad.push_back({std::abs(dJ_dPhi[i]), i});
    }
    std::sort(sorted_grad.rbegin(), sorted_grad.rend());

    // On teste les n_top_cells premières + les n_top_cells du milieu
    std::vector<size_t> cells_to_test;
    for (size_t t = 0; t < std::min(n_top_cells, sorted_grad.size()); t++)
        cells_to_test.push_back(sorted_grad[t].second);
    // Quelques cellules du milieu de la distribution (gradient modéré)
    if (sorted_grad.size() > 2 * n_top_cells) {
        size_t mid = sorted_grad.size() / 2;
        for (size_t t = mid; t < std::min(mid + n_top_cells, sorted_grad.size()); t++)
            cells_to_test.push_back(sorted_grad[t].second);
    }

    if (cells_to_test.empty()) {
        std::cout << "[Taylor] Aucune cellule avec gradient non-nul — test annulé.\n";
        return;
    }

    // ------------------------------------------------------------------
    // 2. J_base sur la trajectoire fixe (référence cohérente avec Enzyme)
    // ------------------------------------------------------------------
    std::vector<double> phi_backup(phi_base, phi_base + N);

    double J_base = this->replay_trajectory_cost(phi_base);
    std::cout << std::scientific << std::setprecision(8);
    std::cout << "J_traj(φ_base) = " << J_base << "\n";

    // Vérification de cohérence avec execute_E
    this->reset_to_start();
    double J_execute = this->execute_E(phi_base);
    std::cout << "J_execute(φ_base) = " << J_execute
              << "  (doit être proche de J_traj si la trajectoire est convergée)\n\n";

    // ------------------------------------------------------------------
    // 3. CSV de sortie
    // ------------------------------------------------------------------
    std::ofstream csv("taylor_test_checkpointed.csv");
    csv << "cell_id,grad_enzyme,epsilon,J_traj_pert,fd_grad,ratio\n";
    csv << std::scientific << std::setprecision(10);

    // ------------------------------------------------------------------
    // 4. Test par cellule — range d'epsilon fixe [1e-2 .. 1e-8]
    //    (indépendant de |grad| : on veut voir la convergence, pas forcer
    //     un ε adaptatif qui cache les problèmes)
    // ------------------------------------------------------------------
    const int N_EPS = 7;
    const double eps_start = 1e-2;   // ε max = 1% de perturbation de phi

    for (size_t ci = 0; ci < cells_to_test.size(); ci++)
    {
        size_t cell_i  = cells_to_test[ci];
        double grad_i  = dJ_dPhi[cell_i];
        double phi_i0  = phi_backup[cell_i];

        std::cout << "--- Cellule " << cell_i
                  << "  |  ∂J/∂φ[" << cell_i << "] = " << grad_i;
        if (ci < n_top_cells) std::cout << "  (top " << (ci+1) << " par |grad|)";
        else                  std::cout << "  (milieu de distribution)";
        std::cout << "\n    φ_base[" << cell_i << "] = " << phi_i0 << "\n";

        printf("    %-12s  %-20s  %-20s  %-10s\n",
               "epsilon", "J_traj(φ+ε·eᵢ)", "FD = ΔJ/ε", "Ratio");
        std::cout << "    " << std::string(66, '-') << "\n";

        // -------------------------------------------------------------------
        // FIX boundary-cell eps direction:
        // For phi_base = 1.0 (upper boundary), positive eps pushes phi to 1.01+
        // (outside the valid [0,1] domain), and J may have a kink there.
        // Enzyme computes the interior (left) derivative; FD with positive eps
        // measures the exterior (right) derivative — they don't match.
        // Solution: use eps_sign = -1 for cells at the upper boundary, and
        // +1 for cells at the lower boundary. Always perturb *into* the domain.
        // -------------------------------------------------------------------
        double eps_sign = 1.0;
        if (phi_i0 >= 1.0 - 1e-8) {
            eps_sign = -1.0;
            std::cout << "    [INFO] φ_base=" << phi_i0
                      << " (upper boundary) → perturbation négative (eps_sign=-1)\n";
        } else if (phi_i0 <= 1e-8) {
            eps_sign = +1.0;
            std::cout << "    [INFO] φ_base=" << phi_i0
                      << " (lower boundary) → perturbation positive (eps_sign=+1)\n";
        }

        double eps = eps_start;
        bool converged = false;
        double prev_ratio = 0.0;

        for (int k = 0; k < N_EPS; k++)
        {
            phi_base[cell_i] = phi_i0 + eps_sign * eps;

            double J_pert  = this->replay_trajectory_cost(phi_base);
            double fd_grad = (J_pert - J_base) / (eps_sign * eps);
            double ratio   = (std::abs(grad_i) > 1e-40) ? fd_grad / grad_i : 0.0;

            // Détection de la convergence → 1
            if (k > 0 && std::abs(ratio - 1.0) < 0.05) converged = true;

            printf("    %12.4e  %20.12e  %20.12e  %10.6f%s\n",
                   eps, J_pert, fd_grad, ratio,
                   (std::abs(ratio - 1.0) < 0.05 ? "  ✓" : ""));

            csv << cell_i << "," << grad_i << "," << eps << ","
                << J_pert << "," << fd_grad << "," << ratio << "\n";

            phi_base[cell_i] = phi_i0;   // restaurer avant le prochain ε
            prev_ratio = ratio;
            eps *= 0.1;
        }

        if (converged)
            std::cout << "    → GRADIENT CORRECT (ratio → 1 sur au moins 1 valeur de ε)\n";
        else
            std::cout << "    → GRADIENT SUSPECT (ratio ne converge pas vers 1)\n";
        std::cout << "\n";
    }

    // ------------------------------------------------------------------
    // 6. Test directionnel global (d = ∇J / ||∇J||₂)
    //    [J_traj(φ+ε·d) - J_traj(φ)] / (ε · ||∇J||₂)  →  1
    //
    //    Ce test vérifie TOUS les composantes du gradient simultanément.
    //    Si les cellules individuelles sont OK mais le global ne l'est pas,
    //    il y a un problème de couplage entre cellules.
    // ------------------------------------------------------------------
    double norm_grad = 0.0;
    for (size_t i = 0; i < N; i++) norm_grad += dJ_dPhi[i] * dJ_dPhi[i];
    norm_grad = std::sqrt(norm_grad);

    // ε_départ : la cellule avec |grad| max détermine la perturbation max
    double max_abs_grad = sorted_grad.empty() ? 1.0 : sorted_grad[0].first;
    // On veut eps * (max_abs_grad / norm_grad) ≤ 0.01  (1% de perturbation phi)
    double eps_dir_start = 0.01 * norm_grad / max_abs_grad;

    std::cout << "--- Test directionnel global (d = ∇J / ||∇J||₂) ---\n";
    std::cout << "    ||∇J||₂ = " << norm_grad << "\n";
    printf("    %-12s  %-20s  %-20s  %-10s\n",
           "epsilon", "J_traj(φ+ε·d)", "FD / ||∇J||", "Ratio");
    std::cout << "    " << std::string(66, '-') << "\n";

    double eps_dir = eps_dir_start;
    for (int k = 0; k < N_EPS; k++)
    {
        // -------------------------------------------------------------------
        // FIX directional test denominator:
        // When phi is at the boundary (0 or 1), the clamped perturbation may
        // differ from eps_dir * dJ_dPhi[i]/norm_grad.  For phi=1 and grad>0,
        // the clamp zeroes the perturbation for that cell.
        // The correct expected FD slope is  ∇J · actual_delta  (not eps*||∇J||²)
        // so we compute the actual_delta before replaying.
        // -------------------------------------------------------------------
        double dot_grad_actual = 0.0;
        for (size_t i = 0; i < N; i++) {
            double phi_new = phi_backup[i] + eps_dir * (dJ_dPhi[i] / norm_grad);
            if (phi_new < 0.0) phi_new = 0.0;
            if (phi_new > 1.0) phi_new = 1.0;
            double actual_delta = phi_new - phi_backup[i];
            phi_base[i] = phi_new;
            dot_grad_actual += dJ_dPhi[i] * actual_delta;
        }

        double J_pert = this->replay_trajectory_cost(phi_base);
        // ratio = (J_pert - J_base) / (∇J · actual_delta) → 1
        // (correct even when boundary cells are clamped)
        double delta_J = J_pert - J_base;
        double ratio  = (std::abs(dot_grad_actual) > 1e-30) ? delta_J / dot_grad_actual : 0.0;
        double fd     = (eps_dir > 1e-300) ? delta_J / eps_dir : 0.0;  // for display

        printf("    %12.4e  %20.12e  %20.12e  %10.6f%s\n",
               eps_dir, J_pert, fd, ratio,
               (std::abs(ratio - 1.0) < 0.05 ? "  ✓" : ""));

        csv << "global," << norm_grad << "," << eps_dir << ","
            << J_pert << "," << fd << "," << ratio << "\n";

        for (size_t i = 0; i < N; i++) phi_base[i] = phi_backup[i];
        eps_dir *= 0.1;
    }

    csv.close();
    std::cout << "\n[Taylor] Résultats exportés dans 'taylor_test_checkpointed.csv'\n"
              << "========================================\n\n";

    // Restauration finale garantie
    for (size_t i = 0; i < N; i++) phi_base[i] = phi_backup[i];
}
