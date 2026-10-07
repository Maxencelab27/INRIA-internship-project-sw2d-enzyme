#include "swCell.h"

#include <cmath>
#include <numeric> 
#include <string>

#include "../../common/swLog.h"
#include "swNode.h"


swCell::swCell(size_t i_id_sms, size_t i_id_sw2d, std::vector<swNode*> nodes_ptr)
/**
 * @brief Construct a new swCell::swCell object
 * 
 * @param[in] i_id_sms id (SMS numbering style)
 * @param[in] i_id_sw2d id (SW2D numbering style)
 * @param[in] nodes_ptr list of nodes (swNode *) attached to the cell
 */
{
    id_sms = i_id_sms;
    idx = i_id_sw2d;

    this->nodes.clear();
    size_t nb_nodes = nodes_ptr.size();

    for (size_t p = 0; p < nodes_ptr.size(); ++p)
    {
        // add the node to cell nodes
        nodes.push_back(nodes_ptr[p]);

        // add the cell to node cells
        nodes_ptr[p]->addCell(this);

        // add neighbouring nodes to node
        std::vector<swNode *> neigh_nodes(2);
        if (p == 0)
        {
            neigh_nodes[0] = nodes_ptr[nb_nodes - 1];
            neigh_nodes[1] = nodes_ptr[p + 1];
        }
        else if (p == nb_nodes - 1)
        {
            neigh_nodes[0] = nodes_ptr[p - 1];
            neigh_nodes[1] = nodes_ptr[0];
        }
        else
        {
            neigh_nodes[0] = nodes_ptr[p - 1];
            neigh_nodes[1] = nodes_ptr[p + 1];
        }
        nodes_ptr[p]->addNeighNodes(neigh_nodes);
    }
    computeArea();
    computeGravityCenter();
    computePerimeter();
}

void swCell::addEdgeBoundary(swEdgeBoundary* edge)
/**
 * @brief add a boundary edge attached to the cell
 * @param[in] edge pointer to one of the cell edge
 */
{
    edges_bnd.push_back(edge);
}

void swCell::addEdgeIntern(swEdgeIntern* edge)
/**
 * @brief add an intern edge attached to the cell
 * @param[in] edge pointer to the edge
 */
{
    edges_int.push_back(edge);
}

void swCell::addNeigh(swCell *cell)
/**
 * @brief add a neighbouring cell
 * @param[in] cell pointer to a neighbouring cell
 */
{
    neigh_cells.push_back(cell);
}

bool swCell::computeArea(void)
/**
 * @brief Compute the cell area
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    area = 0;
    size_t nb_nodes = nodes.size();
    for (size_t k=0; k<nb_nodes; k++)
    {
        double cp = nodes[k]->getX() * nodes[(k+1)%nb_nodes]->getY() - nodes[(k+1)%nb_nodes]->getX() * nodes[k]->getY();
        area += cp;
    }
    area = std::fabs(area)/2;

    return EXIT_SUCCESS;
}

bool swCell::computeGravityCenter(void)
/**
 * @brief Compute x, y and z coordinates of the cell gravity center
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    x_g = 0;
    y_g = 0;
    z_b = 0;

    for (swNode* node : nodes)
    {
        x_g += node->getX();
        y_g += node->getY();
        z_b += node->getZ();
    }
    const size_t nb_nodes = nodes.size();
    x_g /= nb_nodes;
    y_g /= nb_nodes;
    z_b /= nb_nodes;

    return EXIT_SUCCESS;
}

bool swCell::computePerimeter(void)
/**
 * @brief compute the cell perimeter
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    const size_t nb_node = nodes.size();

    perimeter = 0;
    for (size_t p=0; p<nb_node; p++)
    {
        double x1 = nodes[p]->getX();
        double y1 = nodes[p]->getY();
        double x2 = nodes[(p+1)%nb_node]->getX();
        double y2 = nodes[(p+1)%nb_node]->getY();

        perimeter += std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
    }

    return EXIT_SUCCESS;
}

bool swCell::correctTopography(std::vector<double> corr_zb, const std::string &param_corr)
/**
 * @brief correct cell topography (zb)
 * 
 * @param[in] corr_zb list of topographies to compute the new one
 * @param[in] param_corr method to apply to compute new topography from the zb list
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    if (param_corr == "min")
    {
        z_b = *std::min_element(corr_zb.begin(), corr_zb.end());
    }
    else if (param_corr == "max") 
    {
        z_b = *std::max_element(corr_zb.begin(), corr_zb.end());
    }
    else if (param_corr == "mean")
    {
        z_b = std::accumulate(corr_zb.begin(), corr_zb.end(), 0.0) / corr_zb.size();
    }
    else    // GCOVR_EXCL_START
    {
        logError() << "  Unknown correction parameter: " << param_corr;
        logError() << "    (should be min, max, or mean)";
        logFatal(__FILE__, __func__, __LINE__);
    }
    return EXIT_SUCCESS;    // GCOVR_EXCL_STOP
} // GCOVR_EXCL_LINE

double swCell::getArea(void)
/**
 * @brief return the cell area
 */
{
    return area;
}

const std::vector<swEdgeBoundary *> &swCell::getEdgesBoundary(void)
/**
 * @brief return the list of the boundary edge pointers
 */
{
    return edges_bnd;
}

const std::vector<swEdgeIntern *> &swCell::getEdgesIntern(void)
/**
 * @brief return the list of the intern edge pointers
 */
{
    return edges_int;
}

size_t swCell::getIdx(void)
/**
 * @brief return the cell id (sw2d numbering style)
 */
{
    return idx;
}

const std::vector<swCell *>& swCell::getNeigh(void)
/**
 * @brief return the list of neighbouring cell pointers
 */
{
    return neigh_cells;
}

const std::vector<swNode *>& swCell::getNodes(void)
/**
 * @brief return the list of node pointers
 */
{
    return nodes;
}

double swCell::getPerimeter(void)
/**
 * @brief return the cell perimeter
 */
{
    return perimeter;
}

size_t swCell::getSMSIdx(void)
/**
 * @brief return the cell id (SMS numbering style)
 */
{
    return id_sms;
}

double swCell::getXG(void)
/**
 * @brief return x-coordinate of the gravity center
 */
{
    return x_g;
}

double swCell::getYG(void)
/**
 * @brief return y-coordinate of the gravity center
 */
{
    return y_g;
}

double swCell::getZBottom(void)
/**
 * @brief return bottom elevation
 */
{
    return z_b;
}

void swCell::setEdgesBnd(const std::vector<swEdgeBoundary *> &edges_bnd_list)
/**
 * @brief sets boundary edges attached to the current cell
 */
{
    edges_bnd = edges_bnd_list;
}

void swCell::setEdgesInt(const std::vector<swEdgeIntern *> &edges_int_list)
/**
 * @brief sets internal edges attached to the current cell
 */
{
    edges_int = edges_int_list;
}

void swCell::setNeigh(const std::vector<swCell *> &neigh_cells_list)
/**
 * @brief sets the list of all cells located next to the current one
 */
{
    neigh_cells = neigh_cells_list;
}

void swCell::setNodes(const std::vector<swNode *> &nodes_list)
/**
 * @brief sets nodes attached to the current cell
 */
{
    nodes = nodes_list;
}

void swCell::setProperties(size_t idx_sw2d, size_t sms_id,
                           double x, double y, double z,
                           double a, double p,
                           size_t nb_n, size_t nb_bnd, size_t nb_int, size_t nb_nei)
/**
 * @brief sets all main attributes for the current cell
 */
{
    idx = idx_sw2d;
    id_sms = sms_id;
    x_g = x;
    y_g = y;
    z_b = z;
    area = a;
    perimeter = p;
    nodes = std::vector<swNode*>(nb_n);
    edges_bnd = std::vector<swEdgeBoundary*>(nb_bnd);
    edges_int = std::vector<swEdgeIntern*>(nb_int);
    neigh_cells = std::vector<swCell*>(nb_nei);
}