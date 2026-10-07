#include "swBoundaryConditionRatingCurve.h"
#include <cmath>
#include "../../../common/param.h"
#include "../../../common/swLog.h"
#include "../../../common/swString.h"

swBoundaryConditionRating::swBoundaryConditionRating(std::string bc_type) : swBoundaryCondition(bc_type)
{
    // do nothing on purpose
}

swBoundaryConditionRating *swBoundaryConditionRating::create(std::string bc_type)
/**
 * @brief instanciator for a swBoundaryConditionRating (or inheritent) object
 * @param[in] bc_type: type of bc declared
 * @return the pointer to the corresponding child object
 */
{
    swBoundaryConditionRating *rating_curve = nullptr;
    if (bc_type == "interp")
    {
        rating_curve = new swBoundaryConditionRatingInterp(bc_type);
    }
    else if (bc_type == "powerlaw")
    {
        rating_curve = new swBoundaryConditionRatingPowerLaw(bc_type);
    }
    else
    {
        logError() << "Creating rating curve";
        logError() << "  Unknown rating curve type:" << bc_type;
        logError() << "  Allowed rating curve types (keyword type) are:";
        logError() << "    - interp   : discharge is interpolated based on the user-defined rating-curve";
        logError() << "    - powerlaw : discharge is computed based on a powerlaw";
        logFatal(__FILE__, __func__, __LINE__);
    }
    return rating_curve;
}

swBoundaryConditionRatingInterp::swBoundaryConditionRatingInterp(std::string type) : swBoundaryConditionRating(type)
{
    // do nothing on purpose
}

void swBoundaryConditionRatingInterp::computeValBC(const double& z)
/**
 * @brief interp the mass flux based on the free surface elevation
 * @param[in] z: value fo the free surface elevation
 */
{

    // compute the value of the mass flux at the boundary
    if (z <= z0)
    {
        bc_val = 0;
    }
    else
    {
        // loop backward
        while ((z < z_inf) && (rank > 0))
        {
            --rank;
            z_sup = z_inf;
            q_sup = q_inf;
            z_inf = z_tab[rank];
            q_inf = q_tab[rank];
        }

        // loop forward
        while ((z > z_sup) && (rank < z_tab.size() - 1))
        {
            ++rank;
            z_inf = z_sup;
            q_inf = q_sup;
            z_sup = z_tab[rank];
            q_sup = q_tab[rank];
        }

        double delta_z = (z - z_inf) / (z_sup - z_inf);
        bc_val = -(q_inf + (q_sup - q_inf) * delta_z);
    }
}

void swBoundaryConditionRatingInterp::initialise(const std::vector<std::string> &data)
{
    std::string line = data[2];
    std::vector<std::string> list = swString::splitStrOnce(swString::trim(line), swString::regexWhiteSpace(), 1);
    if ((list.size() != 2) || (list[0].compare("z0") != 0))
    {
        logError() << "Impossible to create rating curve boundary - 3rd parameter";
        logError() << data;
        logFatal(__FILE__, __func__, __LINE__);
    }
    bool conversion_success = true;
    z0 = swString::stringToDouble(list[1], conversion_success);
    if (conversion_success == false)
    {
        logError() << "Impossible to create rating curve boundary - 3rd parameter";
        logError() << "  Parameter 'z0' has to be number";
        logError() << data;
        logFatal(__FILE__, __func__, __LINE__);
    }

    line = data[3];
    list = swString::splitStr(swString::trim(line), swString::regexWhiteSpace());
    if ((list.size() < 2) || (list[0].compare("h_val") != 0))
    {
        logError() << "Impossible to create rating curve boundary - 4th parameter";
        logError() << data;
        logFatal(__FILE__, __func__, __LINE__);
    }
    for (size_t i = 1; i < list.size() ; ++i)
    {
        bool conversion_success = true;
        double h = swString::stringToDouble(list[i], conversion_success);
        if (conversion_success == false)
        {
            logError() << "Impossible to create rating curve boundary - 4th parameter";
            logError() << "  Parameter 'h_val' has to be a number list";
            logError() << data;
            logFatal(__FILE__, __func__, __LINE__);
        }
        z_tab.push_back(z0 + h);
    }

    line = data[4];
    list = swString::splitStr(swString::trim(line), swString::regexWhiteSpace());
    if ((list.size() < 2) || (list[0].compare("q_val") != 0))
    {
        logError() << "Impossible to create rating curve boundary - 5th parameter";
        logError() << data;
        logFatal(__FILE__, __func__, __LINE__);
    }
    for (size_t i = 1; i < list.size() ; ++i)
    {
        bool conversion_success = true;
        double q = swString::stringToDouble(list[i], conversion_success);
        if (conversion_success == false)
        {
            logError() << "Impossible to create rating curve boundary - 5th parameter";
            logError() << "  Parameter 'q_val' has to be a number list";
            logError() << data;
            logFatal(__FILE__, __func__, __LINE__);
        }
        q_tab.push_back(q);
    }

    if (z_tab.size() != q_tab.size())
    {
        logError() << "Impossible to create rating curve boundary - 4th and 5th parameters";
        logError() << "  Parameters 'h_val' and 'q_val' have to contain number lists with the same number of element";
        logError() << data;
        logFatal(__FILE__, __func__, __LINE__);
    }
}

swBoundaryConditionRatingPowerLaw::swBoundaryConditionRatingPowerLaw(std::string type) : swBoundaryConditionRating(type)
{
    // do nothing on purpose
}

void swBoundaryConditionRatingPowerLaw::computeValBC(const double& z)
/**
 * @brief interpolate the value of the unit-discharge based on the rating curve
 * @param[in] z: free surface elevation in the neighbouring cell
 */
{
    if (z > z0)
    {
        bc_val = -a * std::exp(b * std::log(z - z0));
    }
    else
    {
        bc_val = 0;
    }
}


void swBoundaryConditionRatingPowerLaw::initialise(const std::vector<std::string> &data)
{
    std::string line = data[2];
    std::vector<std::string> list = swString::splitStrOnce(swString::trim(line), swString::regexWhiteSpace(), 1);
    if ((list.size() != 2) || (list[0].compare("z0") != 0))
    {
        logError() << "Impossible to create rating curve boundary - 3rd parameter";
        logError() << data;
        logFatal(__FILE__, __func__, __LINE__);
    }
    bool conversion_success = true;
    z0 = swString::stringToDouble(list[1], conversion_success);
    if (conversion_success == false)
    {
        logError() << "Impossible to create rating curve boundary - 3rd parameter";
        logError() << "  Parameter 'z0' has to be number";
        logError() << data;
        logFatal(__FILE__, __func__, __LINE__);
    }

    line = data[3];
    list = swString::splitStrOnce(swString::trim(line), swString::regexWhiteSpace(), 1);
    if ((list.size() != 2) || (list[0].compare("a") != 0))
    {
        logError() << "Impossible to create rating curve boundary - 4th parameter";
        logError() << data;
        logFatal(__FILE__, __func__, __LINE__);
    }
    conversion_success = true;
    a = swString::stringToDouble(list[1], conversion_success);
    if ((conversion_success == false) || (a < 0))
    {
        logError() << "Impossible to create rating curve boundary - 4th parameter";
        logError() << "  Parameter 'a' has to be a positive number";
        logError() << data;
        logFatal(__FILE__, __func__, __LINE__);
    }

    line = data[4];
    list = swString::splitStrOnce(swString::trim(line), swString::regexWhiteSpace(), 1);
    if ((list.size() != 2) || (list[0].compare("b") != 0))
    {
        logError() << "Impossible to create rating curve boundary - 5th parameter";
        logError() << data;
        logFatal(__FILE__, __func__, __LINE__);
    }
    conversion_success = true;
    b = swString::stringToDouble(list[1], conversion_success);
    if ((conversion_success == false) || (b < 0))
    {
        logError() << "Impossible to create rating curve boundary - 5th parameter";
        logError() << "  Parameter 'b' has to be a positive number";
        logError() << data;
        logFatal(__FILE__, __func__, __LINE__);
    }

    // check the required parameters have been set
    if ((std::fabs(a) < small_real) ||
        (std::fabs(b) < small_real))
    {
        logError() << "Unsufficient parametrization";
        logError() << "User has to provide z0, a and b such as q=a.(z - z0)^b";
        logError() << "absolute value of a and b have to be greater than" << small_real;
        logFatal(__FILE__, __func__, __LINE__);
    }
}
