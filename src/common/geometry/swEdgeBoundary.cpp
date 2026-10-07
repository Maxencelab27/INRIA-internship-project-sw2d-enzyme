#include "swEdgeBoundary.h"

#include <algorithm>
#include <cmath>

#include "swCell.h"
#include "../../common/swLog.h"
#include "swNode.h"

swEdgeBoundary::swEdgeBoundary(size_t id_sw2d, const std::vector<swNode *> edge_nodes,
                               swCell *edge_cell, int code)
/**
 * @brief Construct a new swEdgeBoundary::swEdgeBoundary object
 * 
 * @param[in] id_sw2d id (SW2D numbering style)
 * @param[in] edge_nodes list of nodes (swNode *) attached to the edge
 * @param[in] edge_cell cell (swCell *) attached to the edge
 * @param[in] bc bc number
 */
{
    idx = id_sw2d;
    nodes = edge_nodes;
    cell = edge_cell;
    bc_code = code;

    // link edge to nodes
    for (swNode* node : nodes)
    {
        node->addEdgeBoundary(this);
    }

    // link edge to cells
    cell->addEdgeBoundary(this);

    computeCoord();
    computeLength();
    computeNormal();
    correctOrientation();
    computeDNorm();
}

bool swEdgeBoundary::computeCoord(void)
/**
 * @brief compute the coordinates of the middle point
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    x_m = (nodes[0]->getX() + nodes[1]->getX()) / 2;
    y_m = (nodes[0]->getY() + nodes[1]->getY()) / 2;
    z_m = (nodes[0]->getZ() + nodes[1]->getZ()) / 2;

    return EXIT_SUCCESS;
}

bool swEdgeBoundary::computeDNorm(void)
/**
 * @brief compute boundary-cell distance in the normal direction to the edge
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    d_norm = norm_x * cell->getXG()
           + norm_y * cell->getYG();
    d_norm = std::fabs(d_norm);

    return EXIT_SUCCESS;
}

bool swEdgeBoundary::computeLength(void)
/**
 * @brief compute the edge length
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    length = sqrt(pow(nodes[0]->getX() - nodes[1]->getX(), 2) +
                  pow(nodes[0]->getY() - nodes[1]->getY(), 2));
    
    return EXIT_SUCCESS;
}

bool swEdgeBoundary::computeNormal(void)
/**
 * @brief compute the components of the unit vector normal to the edge
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    norm_x = (nodes[1]->getY() - nodes[0]->getY()) / length;
    norm_y = -(nodes[1]->getX() - nodes[0]->getX()) / length;

    return EXIT_SUCCESS;
}

bool swEdgeBoundary::correctOrientation(void)
/**
 * @brief check the orientation of the normal vector (from edge to cell)
 * @brief make sure the first node is South (second North)
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    if ((norm_x * (cell->getXG() - x_m) +
         norm_y * (cell->getYG() - y_m)) < 0)
    {
        norm_x = -norm_x;
        norm_y = -norm_y;
        std::reverse(nodes.begin(), nodes.end());
    }
    return EXIT_SUCCESS;
}

int swEdgeBoundary::getBCCode(void)
/**
 * @brief return the bc code of the edge
 */
{
    return bc_code;
}

swCell *swEdgeBoundary::getCell(void)
/**
 * @brief return the cell attached to the edge
 */
{
    return cell;
}

double swEdgeBoundary::getDnorm(void)
/**
 * @brief return the distance boundary-cell in the normal direction to the edge
 */
{
    return d_norm;
}

size_t swEdgeBoundary::getIdx(void)
/**
 * @brief return the edge id in sw2d numbering style (continuous and starting from 0)
 */
{
    return idx;
}

double swEdgeBoundary::getLength(void)
/**
 * @brief return the length of the edge
 */
{
    return length;
}

const std::vector<swNode*>& swEdgeBoundary::getNodes(void)
/**
 * @brief return the edge nodes
 */
{
    return nodes;
}

double swEdgeBoundary::getNormX(void)
/**
 * @brief return the x-component of the normal unit vector
 */
{
    return norm_x;
}

double swEdgeBoundary::getNormY(void)
/**
 * @brief return the y-component of the normal unit vector
 */
{
    return norm_y;
}

double swEdgeBoundary::getXM(void)
/**
 * @brief return the x-coordinate of the edge middle point
 */
{
    return x_m;
}

double swEdgeBoundary::getYM(void)
/**
 * @brief return the y-coordinate of the edge middle point
 */
{
    return y_m;
}

double swEdgeBoundary::getZM(void)
/**
 * @brief return the z-coordinate of the edge middle point
 */
{
    return z_m;
}

void swEdgeBoundary::setBCCode(int code)
/**
 * @brief set the BC code of the edge
 * @param[in] id BC code
 */
{
    bc_code = code;
}

void swEdgeBoundary::setCell(swCell* cell_pointer)
/**
 * @brief defines the cell attached to the current boundary edge
 * @param[in] cell_pointer: pointer to the cell object
 */
{
    cell = cell_pointer;
}

void swEdgeBoundary::setProperties(size_t edge_idx, swNode* node_1, swNode* node_2,
                                   double x, double y, double z,
                                   double d, double l, double nx, double ny,
                                   int code)
/**
 * @brief set the main attributes for one boundary edge
 * @param[in] edge_idx: index of the boundary edge
 * @param[in] node_1: node object defining the first end of the edge
 * @param[in] node_2: node object defining the second end of the edge
 * @param[in] x, y: midpoint coordinates (x&y)
 * @param[in] z: midpoint bottom value
 * @param[in] d: normal distance (in the x-direction) between edge and adjoint cell gravity center 
 * @param[in] l: edge length
 * @param[in] nx, ny: edge normal vector
                                   int code
 */
{
    idx = edge_idx;
    nodes.resize(2);
    nodes[0] = node_1;
    nodes[1] = node_2;
    x_m = x;
    y_m = y;
    z_m = z;
    d_norm = d;
    length = l;
    norm_x = nx;
    norm_y = ny;
    bc_code = code;
}
