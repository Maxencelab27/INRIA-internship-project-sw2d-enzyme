#include "swEdgeIntern.h"

#include <algorithm>
#include <cmath>

#include "swCell.h"
#include "swNode.h"

swEdgeIntern::swEdgeIntern(size_t id_sw2d,
                           const std::vector<swNode *> edge_nodes,
                           const std::vector<swCell *> edge_cells)
/**
 * @brief Construct a new swEdgeIntern::swEdgeIntern object
 * 
 * @param[in] id_sw2d id (SW2D numbering style)
 * @param[in] edge_nodes list of nodes (swNode *) attached to the edge
 * @param[in] edge_cells list of cells (swCell *) attached to the edge
 */
{
    idx = id_sw2d;
    nodes = edge_nodes;
    cells = edge_cells;

    // link edge to nodes
    for (swNode* node : nodes)
    {
        node->addEdgeIntern(this);
    }

    // link edge to cells
    for (swCell* cell : cells)
    {
        cell->addEdgeIntern(this);
    }

    computeCoord();
    computeLength();
    computeNormal();

    // link cells to each other
    cells[0]->addNeigh(cells[1]);
    cells[1]->addNeigh(cells[0]);

    correctOrientation();
    computeDNorm();
}

bool swEdgeIntern::computeCoord(void)
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

bool swEdgeIntern::computeDNorm(void)
/**
 * @brief compute cell-cell distance in the normal direction to the edge
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    d_norm = std::abs(norm_x * (cells[1]->getXG() - cells[0]->getXG())
           + norm_y * (cells[1]->getYG() - cells[0]->getYG()));
    return EXIT_SUCCESS;
}

bool swEdgeIntern::computeLength(void)
/**
 * @brief compute the edge length
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    length = sqrt(pow(nodes[0]->getX() - nodes[1]->getX(), 2) +
                  pow(nodes[0]->getY() - nodes[1]->getY(), 2));
    
    return EXIT_SUCCESS;
}

bool swEdgeIntern::computeNormal(void)
/**
 * @brief compute the components of the unit vector normal to the edge
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    norm_x = (nodes[1]->getY() - nodes[0]->getY()) / length;
    norm_y = -(nodes[1]->getX() - nodes[0]->getX()) / length;

    return EXIT_SUCCESS;
}

bool swEdgeIntern::correctOrientation(void)
/**
 * @brief check the orientation of the normal vector (from edge to cell)
 * @brief make sure the first node is South (second North)
 * @return EXIT_SUCCESS/EXIT_FAILURE
 */
{
    if ((norm_x * (cells[1]->getXG() - cells[0]->getXG()) +
         norm_y * (cells[1]->getYG() - cells[0]->getYG())) < 0)
    {
        norm_x = -norm_x;
        norm_y = -norm_y;
        std::reverse(nodes.begin(), nodes.end());
    }

    return EXIT_SUCCESS;
}

const std::vector<swCell*> &swEdgeIntern::getCells(void)
/**
 * @brief return the cells attached to the edge
 */
{
    return cells;
}

const swCell* swEdgeIntern::getCelL()
/**
 * @brief return the pointer to the left-hand cell
 */
{
    return cells[0];
}

const swCell* swEdgeIntern::getCelR(void)
/**
 * @brief return the pointer to the right-hand cell
 */
{
    return cells[1];
}

double swEdgeIntern::getDnorm(void)
/**
 * @brief return the distance cell-cell in the normal direction to the edge
 */
{
    return d_norm;
}

size_t swEdgeIntern::getIdx(void)
/**
 * @brief return the edge id in sw2d numbering style (continuous and starting from 0)
 */
{
    return idx;
}

double swEdgeIntern::getLength(void)
/**
 * @brief return the length of the edge
 */
{
    return length;
}

const std::vector<swNode*>& swEdgeIntern::getNodes(void)
/**
 * @brief return the edge nodes
 */
{
    return nodes;
}

double swEdgeIntern::getNormX(void)
/**
 * @brief return the x-component of the normal unit vector
 */
{
    return norm_x;
}

double swEdgeIntern::getNormY(void)
/**
 * @brief return the y-component of the normal unit vector
 */
{
    return norm_y;
}

double swEdgeIntern::getXM(void)
/**
 * @brief return the x-coordinate of the edge middle point
 */
{
    return x_m;
}

double swEdgeIntern::getYM(void)
/**
 * @brief return the y-coordinate of the edge middle point
 */
{
    return y_m;
}

double swEdgeIntern::getZM(void)
/**
 * @brief return the z-coordinate of the edge middle point
 */
{
    return z_m;
}

void swEdgeIntern::setCells(std::vector<swCell*> &cells_list)
/**
 * @brief sets the list of (2) cells attached to this internal edge
 */
{
    cells = cells_list;
}

void swEdgeIntern::setProperties(size_t edge_idx, swNode* node_1, swNode* node_2,
                                 double x, double y, double z,
                                 double d, double l, double nx, double ny)
/**
 * @brief set the main attributes for one internal edge
 * @param[in] edge_idx: index of the internal edge
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
}
