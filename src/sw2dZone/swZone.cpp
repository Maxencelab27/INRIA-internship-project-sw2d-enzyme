#include "swZone.h"

#include <cmath>

#include "../common/param.h"
#include "../common/geometry/swCell.h"
#include "../common/geometry/swEdgeBoundary.h"
#include "../common/geometry/swEdgeIntern.h"
#include "../common/swLog.h"
#include "../common/geometry/swMesh.h"
#include "../common/geometry/swNode.h"
#include "../common/swString.h"

swZone::swZone(swMesh *mesh, std::string name, std::vector<double>& m_x_list, std::vector<double>& m_y_list)
/**
 * @brief Construct a new swZone::swZone object
 * 
 * @param[in] mesh pointer to swMesh object 
 * @param[in] name zone name
 * @param[in] m_x_list list of x-coordinates of zone vertices
 * @param[in] m_y_list list of y-coordinates of zone vertices
 */
{
    x_list = m_x_list;
    y_list = m_y_list;
    x_list.push_back(x_list[0]);
    y_list.push_back(y_list[0]);
    m_mesh = mesh;
    m_name = name;
}

bool swZone::addVariable(std::string var, double value, std::vector<std::string> default_var)
/**
 * @brief add a variable which value is defined in the zone
 * 
 * @param[in] var variable name
 * @param[in] value variable value
 * @param[in] default_var list of the defined default variables (to check variable name is in it)
 * @return bool EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    if (swString::contains(default_var, var) == false)
    {
        logError() << "cannot assign parameter" << var;
        logError() << "no default value have been specified for this parameter";
        return EXIT_FAILURE;
    }

    variables.push_back(var);
    values.push_back(value);

    return EXIT_SUCCESS;
}

bool swZone::assignElements()
/**
 * @brief Assign cell, boundary edges and internal edges to the zone
 * 
 * @return EXIT_SUCCESS/EXIT_FAILURE 
 */
{
    for (swCell *cell : m_mesh->getCellsList())
    {
        double xg = cell->getXG();
        double yg = cell->getYG();

        if (isInside(xg, yg) == true)
        {
            cell_list.push_back(cell);
        }
    }

    for (swEdgeBoundary *edge_bnd : m_mesh->getEdgeBoundaryList())
    {
        double xg = edge_bnd->getXM();
        double yg = edge_bnd->getYM();

        if (isInside(xg, yg) == true)
        {
            edge_bnd_list.push_back(edge_bnd);
        }
    }

    for (swEdgeIntern *edge_int : m_mesh->getEdgeInternList())
    {
        double xg = edge_int->getXM();
        double yg = edge_int->getYM();

        if (isInside(xg, yg) == true)
        {
            edge_int_list.push_back(edge_int);
        }
    }

    return EXIT_SUCCESS;
}

swZone *swZone::create(swMesh *mesh, std::string &line)
/**
 * @brief Static methods to create a zone object from a line of zone input file
 * 
 * @param[in] mesh pointer to swMesh object
 * @param[in] line text line from input file containing zone geometric information
 * @return swZone* pointer to a swZone object
 */
{
    swZone *zone = nullptr;

    const std::vector<std::string>& line_split = swString::splitStr(line, swString::regexWhiteSpace());
    std::string zone_nom = line_split[0];
    std::string zone_type = line_split[1];
    std::vector<double> x_list_read;
    std::vector<double> y_list_read;

    if (zone_type == "P")
    {
        // Coordinate specified points
        for (size_t i = 0; i < (line_split.size() - 2) / 2; ++i)
        {
            bool conversion_success = true;

            double x = swString::stringToDouble(line_split[2 * (i + 1)], conversion_success);
            if (conversion_success == false)
            {
                logError() << "line is misformatted";
            }

            double y = swString::stringToDouble(line_split[2 * (i + 1) + 1], conversion_success);
            if (conversion_success == false)
            {
                logError() << "line is misformatted";
            }

            x_list_read.push_back(x);
            y_list_read.push_back(y);
        }
    }
    else if (zone_type == "N")
    {
        // Id specified points
        for (size_t i = 2; i < line_split.size(); ++i)
        {
            bool conversion_success = true;

            size_t node_id = swString::stringToUInt(line_split[i], conversion_success);
            if (conversion_success == false)
            {
                logError() << "line is misformatted";
            }

            //-- Find node
            swNode* node = nullptr;
            for (swNode* n : mesh->getNodesList())
            {
                if (n->getIdx() == node_id)
                {
                    node = n;
                    break;
                }
            };
            if (node == nullptr)
            {
                logError() << "Could not find node for id" << node_id;
                logFatal(__FILE__, __func__, __LINE__);
            }
            //== Find node
            x_list_read.push_back(node->getX());
            y_list_read.push_back(node->getY());
        }
    }
    else
    {
        logError() << "Unknown type of zone";
        logInfo() << "  Recognized type:";
        logInfo() << "   - N: points specified by id";
        logInfo() << "   - P: points specified by coordinate";
        return zone;
    }

    zone = new swZone(mesh, zone_nom, x_list_read, y_list_read);

    return zone;
}

std::vector<swCell*> swZone::getCellList()
/**
 * @brief getter for cell_list attribute
 * 
 * @return std::vector<swCell*> list of cells contained in the zone
 */
{
    return cell_list;
}

std::vector<swEdgeBoundary *> swZone::getEdgeBoundaryList()
/**
 * @brief getter for edge_bnd_list attribute
 * 
 * @return std::vector<swEdge*> list of boundary edges contained in the zone
 */
{
    return edge_bnd_list;
}

std::vector<swEdgeIntern *> swZone::getEdgeInternList()
/**
 * @brief getter for edge_int_list attribute
 * 
 * @return std::vector<swEdge*> list of internal edges contained in zone
 */
{
    return edge_int_list;
}

std::string swZone::getName()
/**
 * @brief getter for m_name attribute
 * 
 * @return std::string name of the zone
 */
{
    return m_name;
}

double swZone::getValue(size_t var_id)
/**
 * @brief Get the value assigned to the zone for one of its variables
 * 
 * @param[in] var_id index of the variable in the zone variables list
 * @return double variable value
 */
{
    return values[var_id];
}

std::vector<std::string> swZone::getVariables()
/**
 * @brief getter for the variables attribute
 * 
 * @return std::vector<std::string> variables defined for this zone
 */
{
    return variables;
}

bool swZone::isInside(double x, double y)
/**
 * @brief Check if a point is inside the zone
 * 
 * @param[in] x x-coordinate of the point
 * @param[in] y y-coordinate of the point
 * @return boolean indicating wether the point is inside the zone or not
 */
{
    double theta = 0.;

    for (size_t i = 0; i < x_list.size() - 1; ++i)
    {
        double x1 = x_list[i];
        double y1 = y_list[i];
        double x2 = x_list[i + 1];
        double y2 = y_list[i + 1];
        double prodvec = (y2 - y) * (x1 - x) -
                         (x2 - x) * (y1 - y);
        double prodscal= (x1 - x) * (x2 - x) +
                         (y1 - y) * (y2 - y);
        double norm1 = std::sqrt((x1 - x) * (x1 - x) +
                                 (y1 - y) * (y1 - y));
        double norm2 = std::sqrt((x2 - x) * (x2 - x) +
                                 (y2 - y) * (y2 - y));
        norm1 = std::max(norm1, small_real);
        norm2 = std::max(norm2, small_real);
        double dtheta = prodvec / (norm1 * norm2);
        double cosine = prodscal / (norm1 * norm2);

        if (cosine > 0.)
        {
            dtheta = asin(dtheta);
        }
        else
        {
            if (dtheta < 0.)
            {
                dtheta = -4.  * atan(1.) - asin(dtheta);
            }
            else
            {
                dtheta = 4. * atan(1.) - asin(dtheta);
            }
        }
        theta = theta + dtheta;
    }

    if (std::fabs(theta) > 6.)
    {
        // Total theta variation is non-zero: point is inside polygon
        return true;
    }
    else
    {
        return false;
    }
}
