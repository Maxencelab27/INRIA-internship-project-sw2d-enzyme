#include <cmath>
#include <vector>

#include "./swDdp.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"
#include "../common/swSource.h"
#include "../../../common/param.h"


void swDdp::sourceFriction(std::string const& model_friction, double const& dt)
/**
 * @brief advances the solution in time for the friction step (Chezy's model)
 * @param[in] model_friction type of friction model (Chezy, Manning, Strickler)
 * @param[in] dt Computational time step
 * @return nothing
 */
{
    if (model_friction == "chezy")
    {
        sourceFrictionChezy(dt);
    }
    else if (model_friction == "manning")
    {
        sourceFrictionManning(dt);
    }
    else if (model_friction == "strickler")
    {
        sourceFrictionStrickler(dt);
    }
    else
    {
        logError() << "Friction model unknown : " << model_friction;
        logFatal(__FILE__, __func__, __LINE__);
    }
}


/** @todo: adapt code */
void swDdp::sourceFrictionChezy(double const& dt)
/**
 * @brief advances the solution in time for the friction step (Chezy's model)
 * @param[in] dt Computational time step
 * @return nothing
 * @par Modified attributes:
 *   - thetaW_velocity_x
 *   - thetaW_velocity_y
 */
{
    ///-- Time splitting discretisation using analytical solution of linearized equation
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        double theta = thetaW[idx_cell];
        if (theta > std::max(h_min, small_real)) // make sure the cell is not dry
        {
            double u = thetaW_velocity_x[idx_cell] / theta;
            double v = thetaW_velocity_y[idx_cell] / theta;   
            double cos_alpha = friction_cos_alpha[idx_cell]; 
            double sin_alpha = friction_sin_alpha[idx_cell];
            // Rotate discharge vector in the friction_alpha direction
            double u_local =  u * cos_alpha + v * sin_alpha; 
            double v_local = -u * sin_alpha + v * cos_alpha;
            double u_norm = std::sqrt(u * u + v * v);
            double chezy_1 = friction_coeff_direction_1[idx_cell];
            double chezy_2 = friction_coeff_direction_2[idx_cell];
            // Apply Chezy's formulae in the principal directions
            double factor_1 = swSource::frictionChezy1D(chezy_1, theta, u_norm, dt);
            double factor_2 = swSource::frictionChezy1D(chezy_2, theta, u_norm, dt);
            u_local *= factor_1;
            v_local *= factor_2;
            // Rotate back in the original frame
            thetaW_velocity_x[idx_cell] = (u_local * cos_alpha - v_local * sin_alpha) * theta;
            thetaW_velocity_y[idx_cell] = (u_local * sin_alpha + v_local * cos_alpha) * theta;
        } 
        else
        {
            thetaW_velocity_x[idx_cell] = 0.;
            thetaW_velocity_y[idx_cell] = 0.;
        }
    }
    //==  Time splitting discretisation using analytical solution
}


/** @todo: adapt code */
void swDdp::sourceFrictionManning(double const& dt)
/**
 * @brief advances the solution in time for the friction step (Chezy's model)
 * @param[in] dt Computational time step
 * @return nothing
 * @par Modified attributes:
 *   - thetaW_velocity_x
 *   - thetaW_velocity_y
 */
{
    ///-- Time splitting discretisation using analytical solution of linearized equation
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        double theta = thetaW[idx_cell];
        if (theta > std::max(h_min, small_real)) // make sure the cell is not dry
        {
            double u = thetaW_velocity_x[idx_cell] / theta;
            double v = thetaW_velocity_y[idx_cell] / theta;
            double cos_alpha = friction_cos_alpha[idx_cell]; 
            double sin_alpha = friction_sin_alpha[idx_cell];
            // Rotate discharge vector in the friction_alpha direction
            double u_local = u * cos_alpha + v * sin_alpha;
            double v_local = -u * sin_alpha + v * cos_alpha;
            double u_norm = std::sqrt(u * u + v * v);
            double manning_1 = friction_coeff_direction_1[idx_cell];
            double manning_2 = friction_coeff_direction_2[idx_cell];
            // Apply Manning's formulae in the principal directions
            double factor_1 = swSource::frictionManning1D(manning_1, theta, u_norm, dt);
            double factor_2 = swSource::frictionManning1D(manning_2, theta, u_norm, dt);
            u_local *= factor_1;
            v_local *= factor_2;
            // Rotate back in the original frame
            thetaW_velocity_x[idx_cell] = (u_local * cos_alpha - v_local * sin_alpha) * theta;
            thetaW_velocity_y[idx_cell] = (u_local * sin_alpha + v_local * cos_alpha) * theta;
        }
        else
        {
            thetaW_velocity_x[idx_cell] = 0.;
            thetaW_velocity_y[idx_cell] = 0.;
        }
    }
    //==  Time splitting discretisation using analytical solution
}


/** @todo: adapt code */
void swDdp::sourceFrictionStrickler(double const& dt)
/**
 * @brief advances the solution in time for the friction step (Chezy's model)
 * @param[in] dt Computational time step
 * @return nothing
 * @par Modified attributes:
 *   - thetaW_velocity_x
 *   - thetaW_velocity_y
 */
{
    ///-- Time splitting discretisation using analytical solution of linearized equation
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        double theta = thetaW[idx_cell];
        if (theta > std::max(h_min, small_real)) // make sure the cell is not dry
        {
            double u = thetaW_velocity_x[idx_cell] / theta;
            double v = thetaW_velocity_y[idx_cell] / theta;
            double cos_alpha = friction_cos_alpha[idx_cell]; 
            double sin_alpha = friction_sin_alpha[idx_cell];
            // Rotate discharge vector in the friction_alpha direction
            double u_local = u * cos_alpha + v * sin_alpha;
            double v_local = -u * sin_alpha + v * cos_alpha;
            double u_norm = std::sqrt(u * u + v * v);
            double strickler_1 = friction_coeff_direction_1[idx_cell];
            double strickler_2 = friction_coeff_direction_2[idx_cell];
            // Apply Strickler's formulae in the principal directions
            double factor_1 = swSource::frictionStrickler1D(strickler_1, theta, u_norm, dt);
            double factor_2 = swSource::frictionStrickler1D(strickler_2, theta, u_norm, dt);
            u_local *= factor_1;
            v_local *= factor_2;
            // Rotate back in the original frame
            thetaW_velocity_x[idx_cell] = (u_local * cos_alpha - v_local * sin_alpha) * theta;
            thetaW_velocity_y[idx_cell] = (u_local * sin_alpha + v_local * cos_alpha) * theta;
        }
        else
        {
            thetaW_velocity_x[idx_cell] = 0.;
            thetaW_velocity_y[idx_cell] = 0.;
        }
    }
    //==  Time splitting discretisation using analytical solution
}


/** @todo: adapt code */
void swDdp::sourceSingularHeadLoss(const double& dt)
/**
 * @brief compute singular head loss source term
 * @param[in] dt Computational time step
 * @return nothing
 * @par Modified attributes:
 *   - thetaW_velocity_x
 *   - thetaW_velocity_y
 */
{
    for (size_t idx_cell = 0; idx_cell < cell_nb; idx_cell++)
    {
        double alpha = singular_head_loss_alpha[idx_cell];
        double cos_alpha = cos(alpha);
        double sin_alpha = sin(alpha);
        double qx =  thetaW_velocity_x[idx_cell] * cos_alpha
                   + thetaW_velocity_y[idx_cell] * sin_alpha;
        double qy = - thetaW_velocity_x[idx_cell] * sin_alpha
                    + thetaW_velocity_y[idx_cell] * cos_alpha;
        double theta = thetaW[idx_cell];
        if (theta > std::max(h_min, small_real))
        {
            double AAx = 0;
            double AAy = 0;
            if (theta > singular_head_loss_h_sing[idx_cell])
            {
                AAx = std::sqrt(qx * qx + qy * qy) * singular_head_loss_coeff_direction_1[idx_cell] / theta;
                AAy = std::sqrt(qx * qx + qy * qy) * singular_head_loss_coeff_direction_2[idx_cell] / theta;
            }
            else
            {
                AAx = 0.;
                AAy = 0.;
            }
            qx = qx * std::exp(-AAx * dt);
            qy = qy * std::exp(-AAy * dt);
        }
        else
        {
            qx = 0.;
            qy = 0.;
        }
        thetaW_velocity_x[idx_cell] = qx * cos_alpha - qy * sin_alpha;
        thetaW_velocity_y[idx_cell] = qx * sin_alpha + qy * cos_alpha;
    }
}