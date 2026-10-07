#pragma once

#include <string>
#include <vector>
#include "../../../common/swBoundaryCondition.h"

class swBoundaryConditionRating : public swBoundaryCondition
/**
 * @brief The swBoundaryConditionRating class
 */
{
public:
    using swBoundaryCondition::swBoundaryCondition;
    swBoundaryConditionRating(std::string type);

    //-- General Methods
    static swBoundaryConditionRating *create(std::string bc_type);
    virtual void initialise(const std::vector<std::string> &data) = 0;
    //== General Methods

protected:
    using swBoundaryCondition::computeValBC;
    virtual void computeValBC(const double&) = 0;

protected:
    //-- Attributes
    std::string bc_type = "q";      /**< type of bc (set to q for reting curve) **/
    double z0 = 0;                  /**< minimum water level below which the discharge is nil**/
    //== Attributes
};

/**
sub-class of swBoundaryConditionRating to handle rating curve described as a series of values
**/
class swBoundaryConditionRatingInterp : public swBoundaryConditionRating
{
public:
    swBoundaryConditionRatingInterp(std::string type);
    ~swBoundaryConditionRatingInterp(void);
    virtual void initialise(const std::vector<std::string> &data) final;

protected:
    //-- General Methods
    virtual void computeValBC(const double& z) final;
    //== General Methods

    //-- Attributes
    size_t rank = 0;        /**< current rank in the rating curve**/
    double z_inf = 0.;      /**< previous value of z (for interpolation) **/
    double z_sup = 0.;      /**< next value of z (for interpolation) **/
    double q_inf = 0.;      /**< previous value of q (for interpolation) **/
    double q_sup = 0.;      /**< next value of q (for interpolation) **/

    std::vector<double> z_tab;  /**< list of the z values for interpolation **/
    std::vector<double> q_tab;  /**< list of the q values for interpolation **/
    //== Attributes
};


/**
sub-class of swBoundaryConditionRating to handle rating curve described as a powerlaw
**/
class swBoundaryConditionRatingPowerLaw : public swBoundaryConditionRating
{
public:
    swBoundaryConditionRatingPowerLaw(std::string type);
    ~swBoundaryConditionRatingPowerLaw(void);

    //-- General Methods
    virtual void initialise(const std::vector<std::string> &data) final;
    //== General Methods

protected:
    //-- General Methods
    virtual void computeValBC(const double& z) final;
    //== General Methods


private:
    //-- attributes
    double a = 0;   /**< multiplicator coefficient for the powerlaw **/
    double b = 0;   /**< exponent coefficient for the powerlaw **/
    //== attributes
};
