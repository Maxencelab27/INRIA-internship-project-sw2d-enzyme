/**
 * @file swSp_E.cpp
 * @brief Enzyme/AD-only members of swSp that have no native SW2D counterpart.
 *
 * Added for automatic differentiation; not part of native SW2D.
 * Declarations live in swSp.h (AD / ENZYME block).
 */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "swSpE.h"
#include "../swParent.h"
#include "../common/swHyp.h"
#include "../common/swParab.h"
#include "../common/swSource.h"
#include "../../../common/param.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"
#include "../../../common/swString.h"
#include "../../common/swSolverParam.h"

// pour enzyme :/*

void swSpE::allocate_as_adjoint(swSpE* primal_solver) {
    // --- 0. RÉCUPÉRATION DES TAILLES (Données immuables) ---
    this->cell_nb      = primal_solver->cell_nb;
    this->edge_int_nb  = primal_solver->edge_int_nb;
    this->edge_bnd_nb  = primal_solver->edge_bnd_nb;
    this->bc_type_nb   = primal_solver->bc_type_nb;
    this->node_nb      = primal_solver->node_nb;

    this->cfl_max      = primal_solver->cfl_max;
    this->v_max        = primal_solver->v_max;
    this->fr_max       = primal_solver->fr_max;

    // Porosity transform settings : the shadow must agree with the primal on the branch
    // and the (constant) normalization scale, so the differentiated path matches.
    this->phi_min        = primal_solver->phi_min;
    this->phi_transform  = primal_solver->phi_transform;
    this->phi_norm_lo    = primal_solver->phi_norm_lo;
    this->phi_norm_hi    = primal_solver->phi_norm_hi;
    this->phi_norm_scale = primal_solver->phi_norm_scale;

    // --- 1. ALLOCATION DES VARIABLES ACTIVES (Initialisées à 0.0 pour accueillir le gradient) ---
    // Ces vecteurs stockeront les "adjoints" (dJ/dVar)
    this->water_depth.assign(this->cell_nb, 0.0);
    this->unit_discharge_x.assign(this->cell_nb, 0.0);
    this->unit_discharge_y.assign(this->cell_nb, 0.0);
    this->swept_data.assign(this->cell_nb, 0.0);
    
    // TRÈS IMPORTANT : phi dans l'adjoint est le gradient dJ/dPhi. Doit être à 0.0.
    this->phi.assign(this->cell_nb, 0.0); 
    
    this->friction_coeff_direction_1.assign(this->cell_nb, 0.0);
    this->friction_coeff_direction_2.assign(this->cell_nb, 0.0);
    this->friction_cos_alpha.assign(this->cell_nb, 0.0);
    this->friction_sin_alpha.assign(this->cell_nb, 0.0);

    // Shadows du head loss singulier (parametres INACTIFS = ne derivent pas de phi).
    // Comme la friction ci-dessus : on alloue des ombres (0.0) de la bonne taille pour
    // qu'Enzyme dispose d'une memoire shadow lors de la lecture dans sourceSingularHeadLoss_E ;
    // les VRAIES valeurs sont lues sur le primal pendant la passe forward de l'autodiff.
    this->singular_head_loss_alpha.assign(this->cell_nb, 0.0);
    this->singular_head_loss_h_sing.assign(this->cell_nb, 0.0);
    this->singular_head_loss_coeff_direction_1.assign(this->cell_nb, 0.0);
    this->singular_head_loss_coeff_direction_2.assign(this->cell_nb, 0.0);

    // Assimilation aux sondes (obs_type=="gauges") : shadows/topologie pour que
    // compute_cost_gauges_ptr soit differentiable sur le solveur adjoint.
    //  - obs_type / n_gauges / gauge_cell_idx / gauge_src_col : copies (constants, indices).
    //  - gauge_w : POIDS constants -> shadow (0.0) de la bonne taille, comme la friction ;
    //    sans lui, Enzyme deref adj->gauge_w vide (nullptr) -> segfault. Sans effet SP/SP.
    this->obs_type       = primal_solver->obs_type;
    this->n_gauges       = primal_solver->n_gauges;
    this->gauge_cell_idx = primal_solver->gauge_cell_idx;
    this->gauge_src_col  = primal_solver->gauge_src_col;
    this->gauge_w.assign(this->n_gauges, 0.0);

    this->bc_val_interpolated.assign(this->edge_bnd_nb, 0.0);

    // --- 2. ALLOCATION DES ÉTATS INTERMÉDIAIRES (Shadows à 0.0) ---
    this->F1_int.assign(this->edge_int_nb, 0.0);
    this->F2_int.assign(this->edge_int_nb, 0.0);
    this->F3_int.assign(this->edge_int_nb, 0.0);
    this->s2L_int.assign(this->edge_int_nb, 0.0);
    this->s2R_int.assign(this->edge_int_nb, 0.0);
    this->lambda_max_int.assign(this->edge_int_nb, 0.0); 
    this->z_riemann_int.assign(this->edge_int_nb, {0.0, 0.0});
    this->q_riemann_int.assign(this->edge_int_nb, {0.0, 0.0});
    this->r_riemann_int.assign(this->edge_int_nb, {0.0, 0.0});
    this->h_riemann_int.assign(this->edge_int_nb, {0.0, 0.0});
    this->lambda_minus_int.assign(this->edge_int_nb, {0.0, 0.0});
    this->lambda_plus_int.assign(this->edge_int_nb, {0.0, 0.0});

    this->F1_bnd.assign(this->edge_bnd_nb, 0.0);
    this->F2_bnd.assign(this->edge_bnd_nb, 0.0);
    this->F3_bnd.assign(this->edge_bnd_nb, 0.0);
    this->s2R_bnd.assign(this->edge_bnd_nb, 0.0);
    this->z_riemann_bnd.assign(this->edge_bnd_nb, 0.0);
    this->q_riemann_bnd.assign(this->edge_bnd_nb, 0.0);
    this->r_riemann_bnd.assign(this->edge_bnd_nb, 0.0);
    this->lambda_max_bnd.assign(this->edge_bnd_nb, 0.0);
    this->lambda_minus_bnd.assign(this->edge_bnd_nb, 0.0);
    this->lambda_plus_bnd.assign(this->edge_bnd_nb, 0.0);
    this->h_riemann_bnd.assign(this->edge_bnd_nb, 0.0); 

    // --- 3. COPIE DE LA TOPOLOGIE (Données de lecture seule - Copie autorisée) ---
    this->cell_area          = primal_solver->cell_area;
    this->cell_z_b           = primal_solver->cell_z_b;
    this->cell_id            = primal_solver->cell_id;
    this->edge_int_celL_idx  = primal_solver->edge_int_celL_idx;
    this->edge_int_celR_idx  = primal_solver->edge_int_celR_idx;
    this->edge_bnd_cell_idx  = primal_solver->edge_bnd_cell_idx;
    this->edge_int_length    = primal_solver->edge_int_length;
    this->edge_bnd_length    = primal_solver->edge_bnd_length;
    this->edge_int_norm_x    = primal_solver->edge_int_norm_x;
    this->edge_int_norm_y    = primal_solver->edge_int_norm_y;
    this->edge_bnd_norm_x    = primal_solver->edge_bnd_norm_x;
    this->edge_bnd_norm_y    = primal_solver->edge_bnd_norm_y;
    this->edge_int_idx       = primal_solver->edge_int_idx;
    this->edge_bnd_idx       = primal_solver->edge_bnd_idx;
    this->edge_total_nb      = primal_solver->edge_total_nb;

    this->cell_edges_int_idx = primal_solver->cell_edges_int_idx;
    this->cell_edges_bnd_idx = primal_solver->cell_edges_bnd_idx;

    this->boussinesq         = primal_solver->boussinesq;
    this->h_min              = primal_solver->h_min;

    // Ponderations du cout (constantes, lues par compute_cost_frame_ptr DANS la tape) :
    // le shadow doit avoir des cartes VALIDES et de meme taille, sinon Enzyme accede a
    // un vecteur vide (taille 0) en differenciant la lecture -> crash hors-bornes.
    this->h_weight_map       = primal_solver->h_weight_map;
    this->q_weight_map       = primal_solver->q_weight_map;

    // FIX C : les facteurs d'arête ne sont PLUS de simples constantes geometriques.
    // computeEdgeFactors_E() les recalcule depuis phi DANS la tape -> ils deviennent
    // des intermediaires ACTIFS. Le shadow (adj) doit donc etre un accumulateur
    // d'adjoint initialise a ZERO (et correctement dimensionne), pas une copie des
    // valeurs du primal (qui polluerait le gradient).
    this->edge_int_factor_1_F1_L.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_1_F1_R.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_2_F2_L.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_2_F2_R.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_2_S2_L.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_2_S2_R.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_2_F3_L.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_2_F3_R.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_3_F2_L.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_3_F2_R.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_3_S2_L.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_3_S2_R.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_3_F3_L.assign(this->edge_int_nb, 0.0);
    this->edge_int_factor_3_F3_R.assign(this->edge_int_nb, 0.0);

    this->edge_bnd_factor_1_F1.assign(this->edge_bnd_nb, 0.0);
    this->edge_bnd_factor_2_F2.assign(this->edge_bnd_nb, 0.0);
    this->edge_bnd_factor_2_S2.assign(this->edge_bnd_nb, 0.0);
    this->edge_bnd_factor_2_F3.assign(this->edge_bnd_nb, 0.0);
    this->edge_bnd_factor_3_F2.assign(this->edge_bnd_nb, 0.0);
    this->edge_bnd_factor_3_S2.assign(this->edge_bnd_nb, 0.0);
    this->edge_bnd_factor_3_F3.assign(this->edge_bnd_nb, 0.0);

    // --- 4. POINTEURS ET TOPOLOGIE FIXE ---
    this->bc_timeseries_list  = primal_solver->bc_timeseries_list;
    this->edge_bnd_bc_pointer = primal_solver->edge_bnd_bc_pointer;
    this->m_io                = primal_solver->m_io;

    this->cell_centre_x       = primal_solver->cell_centre_x;
    this->cell_centre_y       = primal_solver->cell_centre_y;
    this->edge_int_midpoint_x = primal_solver->edge_int_midpoint_x;
    this->edge_int_midpoint_y = primal_solver->edge_int_midpoint_y;
    this->edge_bnd_midpoint_x = primal_solver->edge_bnd_midpoint_x;
    this->edge_bnd_midpoint_y = primal_solver->edge_bnd_midpoint_y;
    this->cell_neigh_cell_idx = primal_solver->cell_neigh_cell_idx;
    this->phi_ptr             = this->phi.data();
    this->current_obs_times   = nullptr;
    this->current_obs_frames  = nullptr;
    this->current_n_obs       = 0;
    
    this->obs_times_vec.clear();
    this->obs_frames_storage.clear();
    this->obs_frames_ptrs.clear();
    
    // --- 5. VARIABLES DE RECONSTRUCTION (Shadows à 0.0) ---
    // Ces variables changent à chaque itération, l'adjoint doit partir de 0.
    this->z_slope_x.assign(this->cell_nb, 0.0);
    this->z_slope_y.assign(this->cell_nb, 0.0);
    this->q_slope_x.assign(this->cell_nb, 0.0);
    this->q_slope_y.assign(this->cell_nb, 0.0);
    this->r_slope_x.assign(this->cell_nb, 0.0);
    this->r_slope_y.assign(this->cell_nb, 0.0);

    this->z_cell.assign(this->cell_nb, 0.0);
    this->dz_dx.assign(this->cell_nb, 0.0);
    this->dz_dy.assign(this->cell_nb, 0.0);

    this->h_init.assign(this->cell_nb, 0.0);
    this->qx_init.assign(this->cell_nb, 0.0);
    this->qy_init.assign(this->cell_nb, 0.0);

    // --- 6. PARAMÈTRES TEMPORELS ---
    this->t_0    = primal_solver->t_0;
    this->t_max  = primal_solver->t_max;
    this->t      = primal_solver->t_0;
    this->dt_max = primal_solver->dt_max;
    this->dt     = primal_solver->dt_max;
    this->n_it_div_max = primal_solver->n_it_div_max;

    // ATTENTION : Suppression des copies directes de z_cell, dz_dx, dz_dy ici.
    // L'adjoint n'a pas besoin des valeurs numériques du primal dans ses membres "actifs",
    // sinon Enzyme confond la valeur et sa dérivée.
}

 double swSpE::specificForce_E(double boussinesq, double phi, double h, double u)
   {
       return phi * (boussinesq * h * u * u + 0.5 * gravity * h * h);
   }

void swSpE::update_phi_norm_bounds(const double* m_map)
/**
 * @brief (Normalize mode only) recompute lo/hi = min/max of the optimization variable m
 *        over the FREE cells (phi_norm_free_idx; all cells if empty), and store them in
 *        phi_norm_lo / phi_norm_hi. Call ONCE per gradient evaluation, BEFORE the Enzyme
 *        tape, so these bounds are constants for that iterate (clean affine derivative).
 */
{
    if (phi_transform != PhiTransform::Normalize) return;

    double lo =  std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    if (!phi_norm_free_idx.empty()) {
        for (size_t idx : phi_norm_free_idx) {
            const double v = m_map[idx];
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
    } else {
        for (size_t i = 0; i < cell_nb; i++) {
            const double v = m_map[i];
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
    }
    phi_norm_lo = lo;
    phi_norm_hi = hi;
    // Scale precompute (1/(hi-lo)) AND degenerate guard, both OUTSIDE the tape, so the
    // in-tape transform stays a pure affine (m-lo)*scale (no division -> no Enzyme nan).
    if (hi - lo < 1e-12) {
        phi_norm_lo    = 0.0;   // degenerate (uniform map) -> phi = (m-0)*1 = m
        phi_norm_scale = 1.0;
        logWarn() << "update_phi_norm_bounds : carte (quasi) uniforme [lo=" << lo
                  << ", hi=" << hi << "] -> normalisation desactivee (phi=m) pour cette iteree.";
    } else {
        phi_norm_scale = 1.0 / (hi - lo);
    }
}

void swSpE::resolve_cost_weight(const std::string& spec, std::vector<double>& out, const char* name)
/**
 * @brief Resout une spec de ponderation du cout : soit un FLOAT (poids uniforme), soit un
 *        CHEMIN de carte (CSV cell_id,phi OU SP Unif/Distrib, lu par load_phi_csv) -> carte
 *        par cellule. Defaut/echec -> poids 1 partout.
 */
{
    out.assign(cell_nb, 1.0);
    // Essai 1 : la spec est-elle un scalaire (float) ?
    try {
        size_t pos = 0;
        double v = std::stod(spec, &pos);
        // tout le texte consomme (a un espace pres) -> c'est bien un scalaire
        while (pos < spec.size() && std::isspace((unsigned char)spec[pos])) pos++;
        if (pos == spec.size()) {
            out.assign(cell_nb, v);
            logInfo() << "Ponderation cout '" << name << "' = scalaire " << v;
            return;
        }
    } catch (...) { /* pas un float -> carte */ }

    // Essai 2 : la spec est un chemin de carte (load_phi_csv gere CSV et SP Unif/Distrib).
    const std::string path = m_io->absolutePath(spec);
    if (!load_phi_csv(path, out.data())) {
        logError() << "resolve_cost_weight : echec lecture carte '" << name << "' (" << path
                   << ") -> poids uniforme 1.";
        out.assign(cell_nb, 1.0);
        return;
    }
    logInfo() << "Ponderation cout '" << name << "' = carte " << path << " (" << cell_nb << " cellules)";
}

void swSpE::storeReferenceMap(const double& time)
/**
 * @brief Writes a reference/observation map directly into reference/ResultHydro_p<time>.txt
 *        Same "default" column layout as the output hydro maps (x, y, zb, h, qx, qy, |q|),
 *        which is exactly what load_single_observation_file() re-reads (cols 0..5).
 *        Called from the time loop during a gen-ref forward run (gen_ref_mode == true).
 * @param[in] time : current simulation time (must land on a dt_map_ad multiple)
 */
{
    ///-- Output file name : <state_dump_dir>/ResultHydro_p<time>.txt
    //    state_dump_dir = "reference" (gen-ref) ou "output" (mode forward).
    const std::string sep = swString::dirSeparator();
    std::string file_name = state_dump_dir + sep + "ResultHydro_p";
    m_io->computeFileNameFromTime(time, file_name);   // append time + ".txt", absolutise
    //== Reference file name

    ///-- Build the same columns as map_storage_mode == "default"
    std::vector<std::string> param_name =
        { "x", "y", "zb", "water_depth", "unit_discharge_x", "unit_discharge_y", "unit_discharge_norm" };
    std::vector<std::vector<double>> param_array;
    param_array.resize(param_name.size());
    for (size_t n = 0; n < cell_nb; n++)
    {
        param_array[0].push_back(cell_centre_x[n]);
        param_array[1].push_back(cell_centre_y[n]);
        param_array[2].push_back(cell_z_b[n]);
        param_array[3].push_back(water_depth[n]);
        param_array[4].push_back(unit_discharge_x[n]);
        param_array[5].push_back(unit_discharge_y[n]);
        double q_norm = std::sqrt(unit_discharge_x[n] * unit_discharge_x[n]
                                + unit_discharge_y[n] * unit_discharge_y[n]);
        param_array[6].push_back(q_norm);
    }
    //== Build columns

    if (swIO::createTxtFileMapDecimal(file_name, param_name, param_array, n_decimal) != EXIT_SUCCESS)
    {
        logError() << "storeReferenceMap - echec d'ecriture de " << file_name;
        return;
    }
    n_ref_written++;
}
