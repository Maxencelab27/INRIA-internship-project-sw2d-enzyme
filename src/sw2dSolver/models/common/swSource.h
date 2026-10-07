#pragma once
#include <string>
#include <vector>

class swSource
{
public:

    swSource(void);
    ~swSource(void);

    //-- General methods
    static double frictionChezy1D(const double& chezy_coeff, const double& hydraulic_radius, const double& u_norm, const double& dt);
    static double frictionManning1D(const double& manning_coeff, const double& hydraulic_radius, const double& u_norm, const double& dt);
    static double frictionStrickler1D(const double& strickler_coeff, const double& hydraulic_radius, const double& u_norm, const double& dt);
    static void   windConvertFromPolarToCartesian(const std::vector<double>& norm, const std::vector<double>& angle_degrees,
                                                    std::vector<double>& component_x, std::vector<double>& component_y);
    static void   precipitationInterp(const double& t, const std::vector<double>& precipitation_time_vector,
                                        const std::vector<double>& precipitation_rate_vector,
                                        double& precipitation_rate_interpolated);
    static double runoffScs(const double& dt, const double& i_a, const double& s, const double& p, double& h_a, double& h_s);
    //== General methods

    //-- Attributes
    //---- Friction attributes  
    static std::vector<double> friction_alpha;
    static std::vector<double> friction_coeff_direction_1;
    static std::vector<double> friction_coeff_direction_2;
    static std::vector<double> friction_cos_alpha;
    static std::vector<double> friction_sin_alpha;
    //==== Friction attributes

    //---- SCS runoff attributes
    static std::vector<double> scs_i_a;    /**< [i]: SCS model - Initial accumulation in cell i*/
    static std::vector<double> scs_h_a ;   /**< [i]: SCS model - Water depth in the Ia compartment for cell i */
    static std::vector<double> scs_h_s;    /**< [i]: SCS model - Water depth in the S compartment for cell i */
    static std::vector<double> scs_s;      /**< [i]: SCS model - Storage capacity in cell i*/
    //==== SCS runoff attributes

    //== Attributes
    
};