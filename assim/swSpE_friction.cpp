/**
 * @file swSpSourceFriction_E.cpp
 * @brief Differentiated (Enzyme) counterparts of swSpSourceFriction.cpp.
 *
 * Added for automatic differentiation; not part of native SW2D.
 * Declarations live in swSp.h (AD / ENZYME block).
 */

#include <cmath>
#include <vector>

#include "swSpE.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"
#include "../common/swSource.h"
#include "../../../common/param.h"


void swSpE::sourceFrictionStrickler_E(double dt) 
{
    // On définit des constantes locales pour éviter les appels à std::max ou des branches complexes
    const double h_threshold = (h_min > small_real) ? h_min : small_real;
    const double g = 9.81; // Assure-toi que c'est la valeur utilisée dans ton swSource
    const double eps = 1e-10;

    for (size_t i = 0; i < cell_nb; i++) 
    {
        double h = water_depth[i];

        if (h > h_threshold) 
        {
            // 1. Calcul des vitesses locales
            double h_safe = h + small_real;
            double u = unit_discharge_x[i] / h_safe;
            double v = unit_discharge_y[i] / h_safe;
            
            double cos_a = friction_cos_alpha[i]; 
            double sin_a = friction_sin_alpha[i];

            // Rotation vers le repère local de friction
            double u_local = u * cos_a + v * sin_a;
            double v_local = -u * sin_a + v * cos_a;
            
            // Norme de la vitesse (essentiel pour Strickler)
            double u_norm = std::sqrt(u * u + v * v + 1e-16);

            // 2. Récupération des coefficients K (ceux que tu calibres !)
            //double K1 = friction_coeff_direction_1[i];
            //double K2 = friction_coeff_direction_2[i];

            // K1 , K2 safe :
            double K1 = (friction_coeff_direction_1[i] > eps) ? friction_coeff_direction_1[i] : 20.0;
            double K2 = (friction_coeff_direction_2[i] > eps) ? friction_coeff_direction_2[i] : 20.0;

            double h_43 = h * std::cbrt(h + 1e-16); 

            // Calcul des facteurs pour chaque direction
            // On ajoute un très petit epsilon au dénominateur pour la stabilité du gradient
            double denom1 = 1.0 + (dt * g * u_norm) / (K1 * K1 * h_43 + eps );
            double denom2 = 1.0 + (dt * g * u_norm) / (K2 * K2 * h_43 + eps );

            u_local /= denom1;
            v_local /= denom2;

            // 4. Rotation inverse et mise à jour des débits unitaires
            unit_discharge_x[i] = (u_local * cos_a - v_local * sin_a) * h;
            unit_discharge_y[i] = (u_local * sin_a + v_local * cos_a) * h;
        }
        else 
        {
            // Gestion des zones sèches : important pour éviter les gradients NaN
            unit_discharge_x[i] = 0.0;
            unit_discharge_y[i] = 0.0;
        }
    }
}


void swSpE::sourceSingularHeadLoss_E(double dt)
{
    // Version DIFFERENTIABLE Enzyme de swSpE::sourceSingularHeadLoss (perte de charge
    // singuliere : dissipation de quantite de mouvement q *= exp(-AA*dt) sur les
    // cellules a structure singuliere, ex. "ville" du cas Toce). Meme physique que le
    // forward classique, rendue Enzyme-safe (cf. sourceFrictionStrickler_E) :
    //   - pas de std::max au point d'usage : plancher h dans la division (h_safe),
    //   - sqrt gardee (+1e-16) pour eviter la derivee infinie a |q|=0,
    //   - pas de log ni d'E/S dans la tape.
    // Les parametres (alpha, coeff, h_sing) sont des CONSTANTES (ne derivent pas de phi)
    // -> branches sur des constantes, tape coherente. Applique APRES la friction, comme
    // dans execute() classique (friction puis head loss).
    const double h_threshold = (h_min > small_real) ? h_min : small_real;
    for (size_t i = 0; i < cell_nb; i++)
    {
        double alpha = singular_head_loss_alpha[i];
        double cos_a = std::cos(alpha);
        double sin_a = std::sin(alpha);
        // Rotation dans le repere principal du tenseur de perte de charge.
        double qx =  unit_discharge_x[i] * cos_a + unit_discharge_y[i] * sin_a;
        double qy = -unit_discharge_x[i] * sin_a + unit_discharge_y[i] * cos_a;
        double h  = water_depth[i];
        if (h > h_threshold)
        {
            double AAx = 0.0;
            double AAy = 0.0;
            // Le terme n'agit qu'au-dela d'une hauteur seuil (comme le classique).
            if (h > singular_head_loss_h_sing[i])
            {
                double h_safe = h + small_real;                       // jamais 1/0 dans la tape
                double qnorm  = std::sqrt(qx * qx + qy * qy + 1e-16); // derivee finie a |q|=0
                AAx = qnorm * singular_head_loss_coeff_direction_1[i] / h_safe;
                AAy = qnorm * singular_head_loss_coeff_direction_2[i] / h_safe;
            }
            // Solution analytique du splitting (amortissement exponentiel).
            qx = qx * std::exp(-AAx * dt);
            qy = qy * std::exp(-AAy * dt);
        }
        else
        {
            qx = 0.0;
            qy = 0.0;
        }
        // Rotation inverse vers le repere global.
        unit_discharge_x[i] = qx * cos_a - qy * sin_a;
        unit_discharge_y[i] = qx * sin_a + qy * cos_a;
    }
}
