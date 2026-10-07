#include "swSpE.h"
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <map> // Ajouté pour observations_map

// On récupère les constantes Enzyme définies ailleurs
extern int enzyme_dup;
extern int enzyme_const;
extern int enzyme_out;

// 1. Sauvegarde de l'état initial
/*void swSpE::save_initial_state() {
    h_init  = water_depth;
    qx_init = unit_discharge_x;
    qy_init = unit_discharge_y;


    this->h_init  = this->water_depth;  // Copie profonde vector
    this->qx_init = this->unit_discharge_x;
    this->qy_init = this->unit_discharge_y;

    this->z_cell_init = this->z_cell; 


    this->swept_data.assign(this->cell_nb, 0.0);
}*/

void swSpE::save_initial_state() {
    // Sauvegarde des CI réelles (water_depth contient déjà les bonnes valeurs)
    this->h_init  = this->water_depth;
    this->qx_init = this->unit_discharge_x;
    this->qy_init = this->unit_discharge_y;
    if (this->swept_data.size() != this->cell_nb)
        this->swept_data.assign(this->cell_nb, 0.0);
    else
        std::fill(this->swept_data.begin(), this->swept_data.end(), 0.0);
}

// 2. Reset pour Enzyme
/*void swSpE::reset_to_start() {
    water_depth      = h_init;
    unit_discharge_x = qx_init;
    unit_discharge_y = qy_init;


    this->swept_data.assign(this->cell_nb, 0.0);
    this->bc_val_interpolated.assign(this->bc_type_nb, 0.0);

    this->F1_int.assign(this->edge_int_nb, 0.0);
    this->F2_int.assign(this->edge_int_nb, 0.0);
    this->F3_int.assign(this->edge_int_nb, 0.0);

    // Reset des pentes et variables de travail
    for (size_t i = 0; i < z_slope_x.size(); i++) z_slope_x[i] = 0.0;
    for (size_t i = 0; i < z_slope_y.size(); i++) z_slope_y[i] = 0.0;
    for (size_t i = 0; i < q_slope_x.size(); i++) q_slope_x[i] = 0.0;
    for (size_t i = 0; i < q_slope_y.size(); i++) q_slope_y[i] = 0.0;
    for (size_t i = 0; i < r_slope_x.size(); i++) r_slope_x[i] = 0.0;
    for (size_t i = 0; i < r_slope_y.size(); i++) r_slope_y[i] = 0.0;
    //for (size_t i = 0; i < z_cell.size();    i++) z_cell[i]    = 0.0; // !!!mettre au ci
    //for (size_t i = 0; i < dz_dx.size();     i++) dz_dx[i]     = 0.0;//!!! mettre aux CI
    //for (size_t i = 0; i < dz_dy.size();     i++) dz_dy[i]     = 0.0;// on part du principe que la topo ne chang e pas lors de la simue. 

    t  = t_0;
    dt = dt_max;
}*/

void swSpE::reset_to_start() {
    // 1. Restauration de l'état physique initial
    this->water_depth      = this->h_init;
    this->unit_discharge_x = this->qx_init;
    this->unit_discharge_y = this->qy_init;

    // 2. Remise à zéro stricte des pentes MUSCL (évite l'effet "mémoire")
    std::fill(this->z_slope_x.begin(), this->z_slope_x.end(), 0.0);
    std::fill(this->z_slope_y.begin(), this->z_slope_y.end(), 0.0);
    std::fill(this->q_slope_x.begin(), this->q_slope_x.end(), 0.0);
    std::fill(this->q_slope_y.begin(), this->q_slope_y.end(), 0.0);
    std::fill(this->r_slope_x.begin(), this->r_slope_x.end(), 0.0);
    std::fill(this->r_slope_y.begin(), this->r_slope_y.end(), 0.0);

    // 3. Remise à zéro des variables de reconstruction
    std::fill(this->z_cell.begin(), this->z_cell.end(), 0.0);
    std::fill(this->dz_dx.begin(),  this->dz_dx.end(),  0.0);
    std::fill(this->dz_dy.begin(),  this->dz_dy.end(),  0.0);

    // 4. Nettoyage des accumulateurs de flux
    std::fill(this->F1_int.begin(), this->F1_int.end(), 0.0);
    std::fill(this->F2_int.begin(), this->F2_int.end(), 0.0);
    std::fill(this->F3_int.begin(), this->F3_int.end(), 0.0);
    std::fill(this->F1_bnd.begin(), this->F1_bnd.end(), 0.0);
    std::fill(this->F2_bnd.begin(), this->F2_bnd.end(), 0.0);
    std::fill(this->F3_bnd.begin(), this->F3_bnd.end(), 0.0);
    std::fill(this->s2R_bnd.begin(), this->s2R_bnd.end(), 0.0);
    std::fill(this->s2L_int.begin(), this->s2L_int.end(), 0.0);
    std::fill(this->s2R_int.begin(), this->s2R_int.end(), 0.0);
    
    std::fill(this->swept_data.begin(), this->swept_data.end(), 0.0);
    this->bc_val_interpolated.assign(this->edge_bnd_nb, 0.0);

    // Toujours ancrer phi_ptr sur le vecteur membre (jamais nullptr).
    // Cela évite le use-after-free sur un temporaire ET garantit que get_phi(i)
    // est valide dès le premier appel après reset (avant execute_E qui recopi phi_map).
    this->phi_ptr = this->phi.data();

    // 5. Reset temps
    this->t  = this->t_0;
    this->dt = this->dt_max;
}


// Version finale pour l'optimisation de la porosité
double swSpE::run_pure_simulation(double* phi_map) {

    if (this->cell_nb <= 0) return 1e30;

    this->reset_to_start();

    double total_accumulated_cost = 0.0;
    int count = 0;

    // Injection de la porosité
    for (size_t i = 0; i < cell_nb; i++) {
        this->phi[i] = phi_map[i];
    }

    while (t < t_max) {
        computeBC_E(t);
        hypCompute_E(n_it_div_max, dt);
        sourceFrictionStrickler_E(dt);

        // Comparaison aux observations temporelles
        // Note: observations_map doit être déclaré dans swSp.h
        for (auto const& [t_obs, h_obs_vector] : observations_map) {
            double diff_t = std::abs(t - t_obs);
            if (diff_t < (dt * 0.5)) {
                total_accumulated_cost += compute_cost_frame(h_obs_vector);
                count++;
            }
        }
        t += dt;
    }

    return (count > 0) ? (total_accumulated_cost / (double)count) : 0.0;
}

// Wrapper pour Enzyme
static void enzyme_wrapper(swSpE* solver, double* K, double* J_out) {
    *J_out = solver->run_pure_simulation(*K);
}

extern "C" {
    void __enzyme_autodiff(void*, ...);
}

double get_global_gradient(swSpE* solver, double K) {
    double dK = 0.0;
    double J_val = 0.0;
    double dJ = 1.0; 

    solver->init_primal_memory();
    swSpE* shadow_solver = new swSpE(); 
    shadow_solver->allocate_as_adjoint(solver);

    __enzyme_autodiff(
        (void*)enzyme_wrapper,
        enzyme_dup, solver, shadow_solver,
        enzyme_dup,   &K, &dK,
        enzyme_dup,   &J_val, &dJ
    );

    return dK;
}

void swSpE::test_enzyme_integration(swSpE* solver) {
    std::cout << "\n[ENZYME] Debut du test d'integration..." << std::endl;
    
    double K_test = 35.0; 
    double grad_enzyme = get_global_gradient(solver, K_test);
    
    std::cout << "Gradient_enzyme dJ/dK : " << grad_enzyme << std::endl;
}