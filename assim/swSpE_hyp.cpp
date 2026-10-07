#include <algorithm>
#include <cmath>
#include <vector>

#include "swSpE.h"
#include "../common/swSource.h"
#include "../common/swBoundaryConditionTimeSeries.h"
#include "../../../common/swBoundaryCondition.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"
#include "../../../common/param.h"

/*void swSpE::hypBoundaryFlux_E(const double& bR, const double& phiR, const double& zbR, const double& zR, const double& qR, const double& rR,
                           const double& lambda_plus_R, 
                           std::string const& bc_type, const double& bc_val,
                           double& flux1, double& flux2, double& flux3, double& source_right)

    
    // NOTE : permier teste , on utilise seulement bc_type = "z"
    
{
    //-- recompute hydrodynamic values in the computational cell
    double hR = zR - zbR;
    double q1 = (qR < h_min) ? 0.0 : qR; // sécuritées pour enzyme 
    double r1 = (rR < h_min) ? 0.0 : rR;
    double hR_safe = (hR < h_min) ? 0.0 : hR;
  
    //== recompute hydrodynamic values in the computational cell
    double hR_denominator = (hR_safe > small_real) ? hR_safe : small_real;
    //-- compute characteristic variables
    double cR = std::sqrt(gravity * hR_denominator);
    double uR = q1 / hR_denominator;
    double L1 = (0.0 > lambda_plus_R) ? 0.0 : lambda_plus_R;

    // specificForce ok pour enzyme 
    double flux2R = specificForce_E(bR, phiR, hR, uR);

    /// version sécuriésée a  revoire :

    // --- REPRISE SÉCURISÉE DE PHIR POUR LES FRONTIÈRES ---
    double phiR_val = 1.0; // Valeur par défaut hors du domaine

    // On récupère l'index de la maille de gauche (qui existe obligatoirement au bord)
    int celL_id = edge_bnd_cell_idx[i]; 

    // On applique une condition de continuité : la porosité au bord externe est égale à celle interne
    if (celL_id >= 0 && celL_id < this->cell_nb) {
        phiR_val = this->phi_ptr ? this->phi_ptr[celL_id] : this->phi[celL_id];
    }

    // On crée la référence attendue par le reste du code sur notre variable locale stable
    const double& phiR = phiR_val;


    //== compute characteristic variables

    //-- Routine Core
    // permier teste , on utilise seulement bc_type = "z"
    // partie a commenter /*if (bc_type == "fr")
    {
        // Critical section
        if ((bc_val < 0.) | (bc_val > 1.))
        {
            logWarn() << "incorrect value for 'fr'-type boundary condition:" << bc_val;
            logWarn() << "  value has to be between 0 and 1";
            logFatal(__FILE__, __func__, __LINE__);
        }
        double c = (uR - 2. * cR) / (- bc_val - 2.);
        c = std::max(c, 0.);
        double u = - bc_val * c;
        double h = c * c / gravity;
        flux1 = phiR * u * h;
        flux2 = flux2R + (flux1 - q1 * phiR) * L1;
    }
    else if (bc_type == "h")
    {
        // Prescribed h
        double h_bc = std::max(bc_val, 0.);
        if (hR > h_min)
        {
            flux1 = phiR * (q1 + L1 * (h_bc - hR));
            flux2 = flux2R + (flux1 - q1 * phiR) * L1;
        }
        else
        {
            double c = std::sqrt(gravity * h_bc);
            flux1 = phiR * h_bc * c;
            flux2 = specificForce(bR, phiR, h_bc, c);
        }
    }
    //else if (bc_type == "z") 
    //### bc_type == "z" : on impose la hauteur de la surface libre ###
    //{
        // Prescribed z
    double h_bc = std::max(bc_val - zbR, 0.);
    if (hR > h_min)
    {
        flux1 = phiR * (q1 + L1 * (h_bc - hR)); // masse 
        flux2 = flux2R + (flux1 - q1 * phiR) * L1; // quantitée de mouvement 
    }
    else
    {
        double c = std::sqrt(gravity * h_bc);
        flux1 = phiR * h_bc * c;
        flux2 = specificForce(bR, phiR, h_bc, c);
    }


      // /*
    else if (bc_type == "q")
    {
        // Prescribed q
        flux1 = bc_val;
        double h_critical = hypHCritical(phiR, bc_val);
        if (h_critical >= hR)
        {
            // water depth in the domain is low => discharge entering with the critical depth
            flux1 = phiR * 1.5 * std::sqrt(gravity * h_critical) * h_critical;
        }
        flux2 = flux2R + (flux1 - q1 * phiR) * L1;
    }
    else
    {
        logInfo() << "BC Type " << bc_type << " not handled ";
        logFatal(__FILE__, __func__, __LINE__);
    }
    //== Routine Core
    //== get boundary condition information

    //-- compute flux of tangeantial momentum
    if (flux1 > 0)
    {
        flux3 = 0;
    }
    else
    {
        double vR = r1 / std::max(hR, small_real);
        flux3 = flux1 * vR;
    }
    //== compute flux of tangeantial momentum
}*/

void swSpE::hypBoundaryFlux_E(const double& bR, const double& phiR, const double& zbR, const double& zR, const double& qR, const double& rR,
                             const double& lambda_plus_R, 
                             std::string const& bc_type, const double& bc_val,
                             double& flux1, double& flux2, double& flux3, double& source_right)
/* Formulae are designed for a domain on the right-hand side of the boundary.
 *
 * Enzyme-safe port of swSpE::hypBoundaryFlux. AVANT, cette fonction ignorait bc_type et
 * appliquait TOUJOURS la formule "z" -> le débit entrant "q" (hydrogramme amont) n'était
 * jamais injecté -> le forward _E vidait le domaine. On rétablit ici le branchement sur le
 * type : "q" (débit imposé), "fr" (sortie critique), "h"/"z" (hauteur/cote imposée).
 *
 * Différentiabilité Enzyme :
 *  - bc_type / bc_val sont enzyme_const (préenregistrés) ; la seule entrée active est phiR
 *    (et l'état de la cellule zR/qR/rR). Les comparaisons sur bc_type sont donc des branches
 *    sur des constantes -> pas de souci pour la tape.
 *  - Toutes les divisions et sqrt sont gardées (max(.,small_real)) pour que la dérivée
 *    inverse reste finie (pas de 0*inf). h_critical n'est calculé que si bc_val>0 (sinon il
 *    vaut 0 en CONSTANTE, donc aucune dérivée 1/q singulière ne traverse la tape).
 */
{
    source_right = 0.0;   // vestigial dans les deux chemins ; valeur déterministe (AD-neutre)

    //-- état hydrodynamique côté domaine (comme le classique : on annule le débit si cellule sèche)
    double hR   = zR - zbR;
    bool   dry  = (hR < h_min);
    double q1   = dry ? 0.0 : qR;
    double r1   = dry ? 0.0 : rR;
    double hR_e = dry ? 0.0 : hR;                              // hR "asséché" (= hR du classique)
    double hden = (hR_e > small_real) ? hR_e : small_real;     // dénominateur gardé
    double cR   = std::sqrt(gravity * hden);
    double uR   = q1 / hden;
    double L1   = (lambda_plus_R > 0.0) ? lambda_plus_R : 0.0;
    double flux2R = specificForce_E(bR, phiR, hR_e, uR);

    flux1 = 0.0;
    flux2 = flux2R;

    //-- Routine Core : branchement sur le type de condition aux limites
    if (bc_type == "q")
    {
        // Débit unitaire imposé (hydrogramme amont). bc_val = q imposé.
        flux1 = bc_val;
        // Profondeur critique h_c = ((q/phi)/g)^(1/3) = exp((1/3)*log((q/phi)/g)).
        // exp/log sont dérivables par Enzyme ; on ne l'évalue que pour q>0 (sinon log(0)=-inf).
        if (bc_val > small_real)
        {
            double phi_safe   = (phiR > phi_min) ? phiR : phi_min;
            double h_critical = std::exp((1.0 / 3.0) * std::log((bc_val / phi_safe) / gravity));
            if (h_critical >= hR_e)
            {
                // domaine peu profond => le débit entre à la profondeur critique
                double hc = (h_critical > small_real) ? h_critical : small_real;
                flux1 = phiR * 1.5 * std::sqrt(gravity * hc) * h_critical;
            }
        }
        flux2 = flux2R + (flux1 - q1 * phiR) * L1;
    }
    else if (bc_type == "fr")
    {
        // Sortie en régime critique (Froude). bc_val dans [0,1] ; on borne sans logger (tape).
        double fr = (bc_val < 0.0) ? 0.0 : ((bc_val > 1.0) ? 1.0 : bc_val);
        double c  = (uR - 2.0 * cR) / (-fr - 2.0);   // dénominateur constant <= -2, jamais nul
        c = (c > 0.0) ? c : 0.0;
        double u  = -fr * c;
        double h  = c * c / gravity;
        flux1 = phiR * u * h;
        flux2 = flux2R + (flux1 - q1 * phiR) * L1;
    }
    else if (bc_type == "h" || bc_type == "z")
    {
        // Hauteur ("h") ou cote ("z") imposée : h_bc ne diffère que par le datum du fond.
        double h_bc = (bc_type == "z") ? std::max(bc_val - zbR, 0.0) : std::max(bc_val, 0.0);
        if (hR_e > h_min)
        {
            flux1 = phiR * (q1 + L1 * (h_bc - hR_e));
            flux2 = flux2R + (flux1 - q1 * phiR) * L1;
        }
        else
        {
            double c = std::sqrt(gravity * ((h_bc > small_real) ? h_bc : small_real));
            flux1 = phiR * h_bc * c;
            flux2 = specificForce_E(bR, phiR, h_bc, c);
        }
    }
    // else : type inconnu -> flux nul (pas de logFatal dans la tape Enzyme).

    //-- flux de quantité de mouvement tangentielle
    if (flux1 > 0.0)
    {
        flux3 = 0.0;
    }
    else
    {
        double vR = r1 / ((hR_e > small_real) ? hR_e : small_real);
        flux3 = flux1 * vR;
    }
}



// On remplace les strings par des types numériques simples
void swSpE::hypCompute_E(size_t n_it_div_max, double& dt_to_use)
{
    hypComputeLRStates_E();

    // 1. Calculer le pas de temps maximum autorisé par la CFL physique
    double dt_phys = hypComputeTimeStep_E();

    // 2. CORRECTION : On ne réduit dt_to_use QUE si la CFL l'exige.
    // Si dt_to_use (notre cible) est plus petit, on le respecte pour ne pas dépasser l'observation.
    if (dt_phys < dt_to_use) {
        dt_to_use = dt_phys;
    }

    hypComputeFluxes_E();

    // 3. CORRECTION : On utilise le pas de temps réel et synchronisé pour la mise à jour des volumes !
    hypComputeBalance_E(dt_to_use);
    
    for (size_t n = 0; n < cell_nb; ++n)
    {
        if (water_depth[n] < 0.)
        {
            water_depth[n] = 0.;
            unit_discharge_x[n] = 0.;
            unit_discharge_y[n] = 0.;
        }
    }
}
    




void swSpE::computeEdgeFactors_E()
{
    // ---- Arêtes de bord ----
    for (size_t i = 0; i < edge_bnd_nb; i++)
    {
        double w    = edge_bnd_length[i];
        size_t c    = edge_bnd_cell_idx[i];
        double area = cell_area[c];
        double nx   = edge_bnd_norm_x[i];
        double ny   = edge_bnd_norm_y[i];
        // FIX NaN gradient : on planche phi DANS la division (jamais 1/0 dans la tape).
        // Avant : (phi>phi_min)?(w/area/phi):0 -> Enzyme evalue/differencie w/area/phi meme
        // pour phi=0 (cellules fermees/murs) -> 1/0=inf, d(1/phi)/dphi=-inf -> 0*inf=NaN
        // (125 cellules phi=0 donnaient un gradient NaN -> ||grad||=nan -> Taylor casse).
        double phi_c = (phi[c] > phi_min) ? phi[c] : phi_min;
        double coeff = (phi[c] > phi_min) ? (w / area / phi_c) : 0.0;

        edge_bnd_factor_1_F1[i] = + coeff;
        edge_bnd_factor_2_F2[i] = + coeff * nx;
        edge_bnd_factor_2_F3[i] = - coeff * ny;
        edge_bnd_factor_2_S2[i] = + coeff * nx;
        edge_bnd_factor_3_F2[i] = + coeff * ny;
        edge_bnd_factor_3_F3[i] = + coeff * nx;
        edge_bnd_factor_3_S2[i] = + coeff * ny;
    }
    // ---- Arêtes internes ----
    for (size_t i = 0; i < edge_int_nb; i++)
    {
        double w     = edge_int_length[i];
        size_t cL    = edge_int_celL_idx[i];
        size_t cR    = edge_int_celR_idx[i];
        double areaL = cell_area[cL];
        double areaR = cell_area[cR];
        double nx    = edge_int_norm_x[i];
        double ny    = edge_int_norm_y[i];
        // FIX NaN gradient : plancher phi DANS la division (cf. arêtes de bord ci-dessus).
        double phi_cL  = (phi[cL] > phi_min) ? phi[cL] : phi_min;
        double phi_cR  = (phi[cR] > phi_min) ? phi[cR] : phi_min;
        double coeff_L = (phi[cL] > phi_min) ? (w / areaL / phi_cL) : 0.0;
        double coeff_R = (phi[cR] > phi_min) ? (w / areaR / phi_cR) : 0.0;

        edge_int_factor_1_F1_L[i] = - coeff_L;
        edge_int_factor_1_F1_R[i] = + coeff_R;
        edge_int_factor_2_F2_L[i] = - coeff_L * nx;
        edge_int_factor_2_F2_R[i] = + coeff_R * nx;
        edge_int_factor_2_F3_L[i] = + coeff_L * ny;
        edge_int_factor_2_F3_R[i] = - coeff_R * ny;
        edge_int_factor_2_S2_L[i] = + coeff_L * nx;
        edge_int_factor_2_S2_R[i] = + coeff_R * nx;
        edge_int_factor_3_F2_L[i] = - coeff_L * ny;
        edge_int_factor_3_F2_R[i] = + coeff_R * ny;
        edge_int_factor_3_F3_L[i] = - coeff_L * nx;
        edge_int_factor_3_F3_R[i] = + coeff_R * nx;
        edge_int_factor_3_S2_L[i] = + coeff_L * ny;
        edge_int_factor_3_S2_R[i] = + coeff_R * ny;
    }
}


void swSpE::hypComputeBalance_E(double dt) // Passage par valeur plus simple pour Enzyme
{
    // ---- 1. Boundary edges ----
    // On suppose que tous les tableaux edge_bnd_factor_... sont pré-alloués
    for (size_t edge_idx = 0; edge_idx < edge_bnd_nb; edge_idx++)
    {
        size_t id_cell = edge_bnd_cell_idx[edge_idx];
        
        water_depth[id_cell] += edge_bnd_factor_1_F1[edge_idx] * F1_bnd[edge_idx] * dt;
        
        unit_discharge_x[id_cell] += ( edge_bnd_factor_2_F2[edge_idx] * F2_bnd[edge_idx]
                                      + edge_bnd_factor_2_S2[edge_idx] * s2R_bnd[edge_idx]
                                      + edge_bnd_factor_2_F3[edge_idx] * F3_bnd[edge_idx]) * dt;
                                      
        unit_discharge_y[id_cell] += ( edge_bnd_factor_3_F2[edge_idx] * F2_bnd[edge_idx]
                                      + edge_bnd_factor_3_S2[edge_idx] * s2R_bnd[edge_idx]
                                      + edge_bnd_factor_3_F3[edge_idx] * F3_bnd[edge_idx]) * dt;
    }

    // ---- 2. Internal edges ----
    for (size_t edge_idx = 0; edge_idx < edge_int_nb ; edge_idx++)
    {
        size_t id_L = edge_int_celL_idx[edge_idx];
        size_t id_R = edge_int_celR_idx[edge_idx];
        
        // Mise à jour h
        double d_h = F1_int[edge_idx] * dt;
        water_depth[id_L] += edge_int_factor_1_F1_L[edge_idx] * d_h;
        water_depth[id_R] += edge_int_factor_1_F1_R[edge_idx] * d_h;
        
        // Mise à jour qx
        unit_discharge_x[id_L] += ( edge_int_factor_2_F2_L[edge_idx] * F2_int[edge_idx]
                                   + edge_int_factor_2_S2_L[edge_idx] * s2L_int[edge_idx]
                                   + edge_int_factor_2_F3_L[edge_idx] * F3_int[edge_idx]) * dt;
                                   
        unit_discharge_x[id_R] += ( edge_int_factor_2_F2_R[edge_idx] * F2_int[edge_idx]
                                   + edge_int_factor_2_S2_R[edge_idx] * s2R_int[edge_idx]
                                   + edge_int_factor_2_F3_R[edge_idx] * F3_int[edge_idx]) * dt;
                                   
        // Mise à jour qy
        unit_discharge_y[id_L] += ( edge_int_factor_3_F2_L[edge_idx] * F2_int[edge_idx]
                                   + edge_int_factor_3_S2_L[edge_idx] * s2L_int[edge_idx]
                                   + edge_int_factor_3_F3_L[edge_idx] * F3_int[edge_idx]) * dt;
                                   
        unit_discharge_y[id_R] += ( edge_int_factor_3_F2_R[edge_idx] * F2_int[edge_idx]
                                   + edge_int_factor_3_S2_R[edge_idx] * s2R_int[edge_idx]
                                   + edge_int_factor_3_F3_R[edge_idx] * F3_int[edge_idx]) * dt;
    }

    // ---- 3. Nettoyage des erreurs d'arrondi ----
    // Attention : mathématiquement, std::fabs crée un "pli" à zéro. 
    // Pour Enzyme, c'est acceptable, mais restez vigilant si le gradient semble instable près de zéro.
    for (size_t n = 0; n < cell_nb; ++n)
    {
        if (std::abs(water_depth[n]) < small_real)      water_depth[n] = 0.0;
        if (std::abs(unit_discharge_x[n]) < small_real) unit_discharge_x[n] = 0.0;
        if (std::abs(unit_discharge_y[n]) < small_real) unit_discharge_y[n] = 0.0;
    }
}

void swSpE::hypComputeFluxes_E(void)
{
    // NOTE : Assurez-vous que topo_discretization == "discontinuous" dans l'input
    
    // -- 1. Boundary edges --
    // bc->computeValBC(t) supprimé : bc_val_interpolated est rempli avant l'appel (par computeBC_E
    // ou par la table préalculée dans le chemin Enzyme). bc_val n'est pas lu ici.

    // boundary edges 
    //core 
    
    for (int i = 0; i < this->edge_bnd_nb; ++i) 
    {
        int cellL = this->edge_bnd_cell_idx[i];

        // 1. Extraction sécurisée par valeur de phi
        double safe_phiR = this->phi_ptr ? this->phi_ptr[cellL] : this->phi[cellL];

        // 2. EXTRACTION SÉCURISÉE PAR VALEUR DE bc_val (ANTI-SEGFAULT)
        double safe_bc_val = 0.0;
        if (i < (int)this->bc_val_interpolated.size()) {
            safe_bc_val = this->bc_val_interpolated[i];
        } else {
            // Sécurité si le vecteur de conditions aux limites n'est pas synchrone avec le maillage
            safe_bc_val = 0.0; 
        }

        // 3. Appel : MÊMES arguments que le forward classique (hypComputeFluxes).
        // CORRECTION (4 args de bord qui divergeaient du classique et faussaient la BC) :
        //  - zR  : z_riemann_bnd[i] (état de Riemann de bord, mis à jour par les LRStates)
        //          et NON z_cell[cellL] -> z_cell n'est PAS mis à jour par le schéma godunov
        //          (seul le MUSCL le faisait) -> zR périmé -> hR faux -> la BC "q" basculait
        //          TOUJOURS en profondeur critique -> ~3x trop de débit injecté (sur-mouillage).
        //  - qR/rR : débits TOURNÉS sur la normale d'arête (q_riemann_bnd/r_riemann_bnd =
        //          q*nx+r*ny / -q*ny+r*nx), et NON unit_discharge_x/y (repère global) -> sinon
        //          la BC voit un débit normal faux (momentum/flux2 erronés).
        //  - lambda : lambda_plus_bnd[i] (comme le classique), et NON lambda_max_bnd.
        this->hypBoundaryFlux_E(
            this->boussinesq[cellL],
            safe_phiR,
            this->cell_z_b[cellL],
            this->z_riemann_bnd[i],
            this->q_riemann_bnd[i],
            this->r_riemann_bnd[i],
            this->lambda_plus_bnd[i],
            this->edge_bnd_bc_pointer[i]->getBCType(),
            safe_bc_val,
            this->F1_bnd[i],
            this->F2_bnd[i],
            this->F3_bnd[i],
            this->s2R_bnd[i]
        );
    }

    // -- 2. Internal edges (Cœur du calcul de gradient pour phi) --
    for (size_t edge_idx = 0 ; edge_idx < edge_int_nb ; edge_idx++)
    {
        const size_t cL = edge_int_celL_idx[edge_idx];
        const size_t cR = edge_int_celR_idx[edge_idx];
        
        const double zbL = cell_z_b[cL];
        const double zbR = cell_z_b[cR];
        const double zL  = z_riemann_int[edge_idx][0];
        const double zR  = z_riemann_int[edge_idx][1];

        // La garde dure h > h_min est supprimée : hypHllc_E applique maintenant
        // une ramp cubique C¹ interne qui éteint les flux continûment → J lisse en phi.
        // On appelle HLLC dès que l'une des deux cellules a de l'eau (h > 0).
        // Pour h ≤ 0 (impossible physiquement), HLLC retourne 0 via l'early-return.
        {
            double f1 = 0, f2 = 0, f3 = 0, sL = 0, sR = 0;

            hypHllc_E(boussinesq[cL], boussinesq[cR], get_phi(cL), get_phi(cR), zbL, zbR,
                      zL, q_riemann_int[edge_idx][0], r_riemann_int[edge_idx][0],
                      zR, q_riemann_int[edge_idx][1], r_riemann_int[edge_idx][1],
                      lambda_minus_int[edge_idx][0], lambda_minus_int[edge_idx][1],
                      lambda_plus_int[edge_idx][0], lambda_plus_int[edge_idx][1],
                      f1, f2, f3, sL, sR);

            F1_int[edge_idx]  = f1;
            F2_int[edge_idx]  = f2;
            F3_int[edge_idx]  = f3;
            s2L_int[edge_idx] = sL;
            s2R_int[edge_idx] = sR;
        }
    }
}




void swSpE::hypComputeLRStates_E()
{
    // CRUCIAL : utiliser le MÊME schéma que le forward classique de référence
    // (hypComputeLRStates(num_scheme)). input.solv: num_scheme=godunov (1er ordre).
    // Avant, on appelait TOUJOURS le MUSCL-EVR (2e ordre) -> schéma numérique DIFFÉRENT
    // du classique -> solution divergente (sous-mouillage global + overshoot 2e ordre au
    // point bas -> emballement du débit de la cellule fr aval 183 -> dt collapse).
    if (num_scheme == "godunov")
        hypGodunovComputeLRStates_E();
    else
        hypMusclEvrComputeLRStates_E();
}

// =============================================================================
// Port Enzyme-safe de hypGodunovComputeLRStates (1er ordre).
// États de Riemann = valeurs au centre de cellule projetées sur la normale d'arête
// (aucune reconstruction de pente, aucune dépendance en dt). Écrit dans des tableaux
// pré-dimensionnés (init_primal_memory) -> pas d'allocation/log dans la tape.
// Garde cellule sèche u=0 (h<=h_min) pour borner la vitesse d'onde (sinon dt s'effondre).
// =============================================================================
bool swSpE::hypGodunovComputeLRStates_E(void)
{
    // ---- Arêtes de bord ----
    for (size_t edge_idx = 0; edge_idx < edge_bnd_nb; edge_idx++)
    {
        double nx = edge_bnd_norm_x[edge_idx];
        double ny = edge_bnd_norm_y[edge_idx];
        size_t id_cell = edge_bnd_cell_idx[edge_idx];
        double bR = boussinesq[id_cell];
        double q  = unit_discharge_x[id_cell];
        double r  = unit_discharge_y[id_cell];
        h_riemann_bnd[edge_idx] = water_depth[id_cell];
        z_riemann_bnd[edge_idx] = water_depth[id_cell] + cell_z_b[id_cell];
        q_riemann_bnd[edge_idx] =  q * nx + r * ny;
        r_riemann_bnd[edge_idx] = -q * ny + r * nx;
        double hR = h_riemann_bnd[edge_idx];
        double uR = (hR > h_min) ? q_riemann_bnd[edge_idx] / hR : 0.0;
        double cR = std::sqrt(gravity * hR);
        lambda_minus_bnd[edge_idx] = lambdaMinus(bR, cR, uR);
        lambda_plus_bnd[edge_idx]  = lambdaPlus(bR, cR, uR);
        lambda_max_bnd[edge_idx]   = std::max(std::abs(lambda_plus_bnd[edge_idx]),
                                              std::abs(lambda_minus_bnd[edge_idx]));
    }

    // ---- Arêtes internes ----
    for (size_t edge_idx = 0; edge_idx < edge_int_nb; edge_idx++)
    {
        double nx = edge_int_norm_x[edge_idx];
        double ny = edge_int_norm_y[edge_idx];

        // côté gauche
        size_t id_celL = edge_int_celL_idx[edge_idx];
        double bL = boussinesq[id_celL];
        double q  = unit_discharge_x[id_celL];
        double r  = unit_discharge_y[id_celL];
        h_riemann_int[edge_idx][0] = water_depth[id_celL];
        z_riemann_int[edge_idx][0] = water_depth[id_celL] + cell_z_b[id_celL];
        q_riemann_int[edge_idx][0] =  q * nx + r * ny;
        r_riemann_int[edge_idx][0] = -q * ny + r * nx;
        double hL = h_riemann_int[edge_idx][0];
        double uL = (hL > h_min) ? q_riemann_int[edge_idx][0] / hL : 0.0;
        double cL = std::sqrt(gravity * hL);
        lambda_minus_int[edge_idx][0] = lambdaMinus(bL, cL, uL);
        lambda_plus_int[edge_idx][0]  = lambdaPlus(bL, cL, uL);

        // côté droit
        size_t id_celR = edge_int_celR_idx[edge_idx];
        double bR = boussinesq[id_celR];
        q = unit_discharge_x[id_celR];
        r = unit_discharge_y[id_celR];
        h_riemann_int[edge_idx][1] = water_depth[id_celR];
        z_riemann_int[edge_idx][1] = water_depth[id_celR] + cell_z_b[id_celR];
        q_riemann_int[edge_idx][1] =  q * nx + r * ny;
        r_riemann_int[edge_idx][1] = -q * ny + r * nx;
        double hR = h_riemann_int[edge_idx][1];
        double uR = (hR > h_min) ? q_riemann_int[edge_idx][1] / hR : 0.0;
        double cR = std::sqrt(gravity * hR);
        lambda_minus_int[edge_idx][1] = lambdaMinus(bR, cR, uR);
        lambda_plus_int[edge_idx][1]  = lambdaPlus(bR, cR, uR);

        lambda_max_int[edge_idx] = std::max({std::abs(lambda_minus_int[edge_idx][0]),
                                             std::abs(lambda_plus_int[edge_idx][0]),
                                             std::abs(lambda_minus_int[edge_idx][1]),
                                             std::abs(lambda_plus_int[edge_idx][1])});
    }
    return EXIT_SUCCESS;
}




double swSpE::hypComputeTimeStep_E(void) {
    double dt = dt_max;
    size_t n = cell_nb;

    // 1. Récupération des pointeurs bruts (très important pour Enzyme)
    double* h = water_depth.data();
    double* swept = this->swept_data.data();
    double* areas = cell_area.data();

    // 2. Réinitialisation propre
    for (size_t i = 0; i < n; i++) {
        swept_data[i] = 0.0;
    }

    // 2.2 autre méthode de réinitialisation 
    //this->swept_member.assign(n, 0.0); 

    // 3. Arêtes limites
    for (size_t edge_idx = 0; edge_idx < edge_bnd_nb; edge_idx++) {
        size_t id_cell = edge_bnd_cell_idx[edge_idx];
        if (h[id_cell] > h_min) {
            swept[id_cell] += lambda_max_bnd[edge_idx] * edge_bnd_length[edge_idx];
        }
    }

    // 4. Arêtes internes
    for (size_t edge_idx = 0; edge_idx < edge_int_nb; edge_idx++) {
        size_t id_L = edge_int_celL_idx[edge_idx];
        size_t id_R = edge_int_celR_idx[edge_idx];
        double length = edge_int_length[edge_idx];
        double lambda = lambda_max_int[edge_idx];

        if (h[id_L] > h_min) swept[id_L] += lambda * length;
        if (h[id_R] > h_min) swept[id_R] += lambda * length;
    }

    // 5. Calcul du dt (SORTI des boucles précédentes)
    for (size_t i = 0; i < n; i++) {
        size_t id = cell_id[i];
        // std::max est ok pour enzyme
        double den = std::max(swept[id], small_real);
        double cfl_dt = cfl_max * (areas[id] / den);
        
        if (cfl_dt < dt) dt = cfl_dt;
    }
/// fonction qui casse touyt!!!!!!
// cette fonction permet d'arriver parfaitement sur le t_max souhaité sur le dernier pas de temps
    // 6. Synchronisation (évite std::floor si possible, ou utilise la version simple)
    /*if (dt_map > 1e-12) {
        double t_next_map = std::floor(t / dt_map + 1.0 + 1e-9) * dt_map;
        if ((t + dt) > t_next_map) dt = t_next_map - t;
    }*/

    if ((t + dt) > t_max) dt = t_max - t;
    if (dt < 1e-6) dt = 1e-6;

    return dt;
}

/* ============================================================
   ANCIENNE VERSION de hypHllc_E — conservée pour comparaison
   (ne pas décommenter : la nouvelle version ci-dessous remplace)

void swSpE::hypHllc_E_OLD(const double& bL, const double& bR,
                   const double& phiL, const double& phiR,
                   const double& zbL,  const double& zbR,
                   const double& zL,   const double& qL, const double& rL,
                   const double& zR,   const double& qR, const double& rR,
                   const double& lambda_minus_L, const double& lambda_minus_R,
                   const double& lambda_plus_L,  const double& lambda_plus_R,
                   double& flux1, double& flux2, double& flux3,
                   double& source_left, double& source_right)
{
    double hL = (zL - zbL > 0.0) ? (zL - zbL) : 0.0;
    double hR = (zR - zbR > 0.0) ? (zR - zbR) : 0.0;

    bool flow_in_both_cells    = false;
    bool flow_in_lhs_cell_alone = false;
    bool flow_in_rhs_cell_alone = false;
    bool flow_nowhere           = false;

    // Régime "rien ne coule" : les deux cellules sont sèches (h ≤ h_min)
    // ou les deux porosités sont nulles
    if (((hL <= h_min) && (hR <= h_min)) || ((phiL <= phi_min) && (phiR <= phi_min)))
    {
        flow_nowhere = true;
    }
    else
    {
        flow_in_both_cells = true;
        // L'eau ne peut pas entrer dans la cellule droite
        if ((phiR <= phi_min) || ((phiR > phi_min) && (hR <= h_min) && (zL <= zbR + h_min)))
        { flow_in_lhs_cell_alone = true; flow_in_both_cells = false; }
        // L'eau ne peut pas entrer dans la cellule gauche
        if ((phiL <= phi_min) || ((phiL > phi_min) && (hL <= h_min) && (zR <= zbL + h_min)))
        { flow_in_rhs_cell_alone = true; flow_in_both_cells = false; }
    }

    if (flow_in_both_cells)
    {
        double hL_safe = (hL > small_real) ? hL : small_real;
        double hR_safe = (hR > small_real) ? hR : small_real;
        double uL = qL / hL_safe;
        double uR = qR / hR_safe;
        double L1 = (lambda_minus_L < lambda_minus_R) ? lambda_minus_L : lambda_minus_R;
        if (L1 > 0.0) L1 = 0.0;
        double L2 = (lambda_plus_L > lambda_plus_R) ? lambda_plus_L : lambda_plus_R;
        if (L2 < 0.0) L2 = 0.0;
        double demoninator = L2 - L1;
        if (demoninator < small_real) demoninator = small_real;
        double phiLR = (phiL < phiR) ? phiL : phiR;
        flux1 = phiLR * (L2 * qL - L1 * qR + L1 * L2 * (zR - zL)) / demoninator;
        double fact = 0.5 * (phiL < phiR ? phiL : phiR) * (hL + hR);
        double s2 = (phiR * hR * hR - phiL * hL * hL) * 0.5 + fact * (zL - zR);
        s2 = s2 * gravity / demoninator;
        source_left  = -L1 * s2;
        source_right =  L2 * s2;
        double FL = specificForce(bL, phiL, hL, uL);
        double FR = specificForce(bR, phiR, hR, uR);
        flux2 = (L2 * (FL + source_left) - L1 * (FR - source_right)
                 + L1 * L2 * (phiR * qR - phiL * qL)) / std::max(L2 - L1, small_real);
        double vL = rL / hL_safe;
        double vR = rR / hR_safe;
        flux3 = (flux1 >= 0) ? (bL * flux1 * vL) : (bR * flux1 * vR);
    }
    else if (flow_in_rhs_cell_alone)
    {
        // Cellule gauche sèche, niveau de la droite sous zbL → remplissage impossible
        flux1 = 0.; flux2 = 0.; flux3 = 0.;
        source_left  = 0.;
        source_right = phiR * 0.5 * gravity * hR * hR;
    }
    else if (flow_in_lhs_cell_alone)
    {
        // Cellule droite sèche, niveau de la gauche sous zbR → remplissage impossible
        flux1 = 0.; flux2 = 0.; flux3 = 0.;
        source_left  = -phiL * 0.5 * gravity * hL * hL;
        source_right = 0.;
    }
    else // flow_nowhere
    {
        flux1 = 0.; flux2 = 0.; flux3 = 0.;
        source_left = 0.; source_right = 0.;
    }
}
   ============================================================ */

void swSpE::hypHllc_E(const double& bL, const double& bR,
                   const double& phiL, const double& phiR, 
                   const double& zbL, const double& zbR,
                   const double& zL, const double& qL, const double& rL,
                   const double& zR, const double& qR, const double& rR,
                   const double& lambda_minus_L, const double& lambda_minus_R,
                   const double& lambda_plus_L, const double& lambda_plus_R, 
                   double& flux1, double& flux2, double& flux3,
                   double& source_left, double& source_right)

/*Pour Enzyme supression des log et des appels a  std min, std max
supression des long commentaires, voir dans swSphyp.cpp pour les commentaires*/

{
    // PORT EXACT du forward classique swSpE::hypHllc (4 régimes durs : flow_nowhere /
    // flow_in_both_cells / flow_in_lhs_cell_alone / flow_in_rhs_cell_alone).
    //
    // On abandonne TEMPORAIREMENT la rampe C¹ lisse (introduite pour la différentiabilité
    // Enzyme) au profit de la FIDÉLITÉ + STABILITÉ du forward : la rampe annulait la source
    // hydrostatique de PAROI (±phi*0.5*g*h²) que le classique applique sur la cellule
    // MOUILLÉE dans les régimes bloqués (lhs_alone/rhs_alone). Sans cette réflexion du
    // momentum, le débit des points bas aval s'accumulait -> vitesse d'onde -> dt collapse.
    //
    // Enzyme-safe : aucun log/resize, divisions gardées (small_real). Les branches booléennes
    // sur h_min rendent J(phi) NON-C¹ (discontinuités) -> à relisser plus tard pour le
    // gradient ; ici la priorité est un forward cohérent avec le classique.
    double hL = (zL - zbL > 0.0) ? (zL - zbL) : 0.0;
    double hR = (zR - zbR > 0.0) ? (zR - zbR) : 0.0;

    bool flow_in_both_cells     = false;
    bool flow_in_lhs_cell_alone = false;
    bool flow_in_rhs_cell_alone = false;
    bool flow_nowhere           = false;

    if (((hL <= h_min) && (hR <= h_min)) || ((phiL <= phi_min) && (phiR <= phi_min)))
    {
        flow_nowhere = true;
    }
    else
    {
        flow_in_both_cells = true;
        // L'eau ne peut pas entrer dans la cellule DROITE -> écoulement dans la gauche seule
        if ((phiR <= phi_min) || ((phiR > phi_min) && (hR <= h_min) && (zL <= zbR + h_min)))
        { flow_in_lhs_cell_alone = true; flow_in_both_cells = false; }
        // L'eau ne peut pas entrer dans la cellule GAUCHE -> écoulement dans la droite seule
        if ((phiL <= phi_min) || ((phiL > phi_min) && (hL <= h_min) && (zR <= zbL + h_min)))
        { flow_in_rhs_cell_alone = true; flow_in_both_cells = false; }
    }

    if (flow_in_both_cells)
    {
        double hL_safe = (hL > small_real) ? hL : small_real;
        double hR_safe = (hR > small_real) ? hR : small_real;
        double uL = qL / hL_safe;
        double uR = qR / hR_safe;
        double L1 = (lambda_minus_L < lambda_minus_R) ? lambda_minus_L : lambda_minus_R;
        if (L1 > 0.0) L1 = 0.0;
        double L2 = (lambda_plus_L > lambda_plus_R) ? lambda_plus_L : lambda_plus_R;
        if (L2 < 0.0) L2 = 0.0;
        double denominator = L2 - L1;
        if (denominator < small_real) denominator = small_real;
        // phiLR = min(phiL, phiR), mais en version LISSE (C1) pour un sous-gradient
        // SYMETRIQUE sous Enzyme. Le min dur "(phiL<phiR)?phiL:phiR" renvoie phiR quand
        // phiL==phiR (cas de l'ebauche uniforme phi=1 sur TOUTES les aretes internes) :
        // Enzyme met alors 100% de la sensibilite de phiLR sur la cellule R et 0% sur L
        // -> le gradient des cellules dominantes est sous-estime (~13% au test de Taylor,
        // ratio plateau != 1). Le min lisse  0.5*(a+b - sqrt((a-b)^2 + eps^2))  a une
        // derivee 0.5/0.5 au point d'egalite, vaut min(a,b) a O(eps^2) pres aux vraies
        // interfaces de barrage (phiL=0 vs phiR=1) -> J redevient C1 -> Taylor converge.
        const double phi_smooth_eps = 1e-3;     // largeur de lissage (echelle de phi in [0,1])
        const double dphi_LR = phiL - phiR;
        double phiLR = 0.5 * (phiL + phiR
                              - std::sqrt(dphi_LR * dphi_LR + phi_smooth_eps * phi_smooth_eps));
        flux1 = phiLR * (L2 * qL - L1 * qR + L1 * L2 * (zR - zL)) / denominator;
        double fact = 0.5 * phiLR * (hL + hR);
        double s2   = (phiR * hR * hR - phiL * hL * hL) * 0.5 + fact * (zL - zR);
        s2 = s2 * gravity / denominator;
        source_left  = -L1 * s2;
        source_right =  L2 * s2;
        double FL = specificForce_E(bL, phiL, hL, uL);
        double FR = specificForce_E(bR, phiR, hR, uR);
        flux2 = (L2 * (FL + source_left) - L1 * (FR - source_right)
                 + L1 * L2 * (phiR * qR - phiL * qL)) / denominator;
        double vL = rL / hL_safe;
        double vR = rR / hR_safe;
        flux3 = (flux1 >= 0.0) ? (bL * flux1 * vL) : (bR * flux1 * vR);
    }
    else if (flow_in_rhs_cell_alone)
    {
        // Cellule R mouillée, L sèche et bloquée (niveau de R sous zbL) : remplissage
        // impossible -> pression hydrostatique de paroi sur R (réflexion du momentum).
        flux1 = 0.0; flux2 = 0.0; flux3 = 0.0;
        source_left  = 0.0;
        source_right = phiR * 0.5 * gravity * hR * hR;
    }
    else if (flow_in_lhs_cell_alone)
    {
        // Cellule L mouillée, R sèche et bloquée : paroi hydrostatique sur L.
        flux1 = 0.0; flux2 = 0.0; flux3 = 0.0;
        source_left  = -phiL * 0.5 * gravity * hL * hL;
        source_right = 0.0;
    }
    else // flow_nowhere (+ cas "impossible" du classique : tout à 0, sans logFatal dans la tape)
    {
        flux1 = 0.0; flux2 = 0.0; flux3 = 0.0;
        source_left = 0.0; source_right = 0.0;
    }
}
