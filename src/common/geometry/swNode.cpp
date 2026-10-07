#include "swNode.h"

#include "../../common/swLog.h"

swNode::swNode(size_t sms_idx, size_t sw2d_idx,
               double x_read, double y_read, double z_read)
/**
 * @brief constructor
 * @param[in] sms_idx node id in SMS (discontinuous numbering and starting from 1)
 * @param[in] sw2d_idx node id in sw2d (continuous numbering and starting from 0)
 * @param[in] x_read x-coordinate
 * @param[in] y_read y-coordinate
 * @param[in] z_read z-coordinate
 */
{
    idx_sms = sms_idx;
    idx_sw2d = sw2d_idx;
    x = x_read;
    y = y_read;
    z = z_read;
    cells.clear();
    edges_bnd.clear();
    edges_int.clear();
    nodes.clear();
}

void swNode::addCell(swCell* cell)
/**
 * @brief add a cell attached to the node
 * @param[in] cell pointer to a neighbouring cell
 */
{
    cells.push_back(cell);
}

void swNode::addEdgeBoundary(swEdgeBoundary* edge)
/**
 * @brief add a boundary edge attached to the node
 * @param[in] edge pointer to one of the cell edge
 */
{
    edges_bnd.push_back(edge);
}

void swNode::addEdgeIntern(swEdgeIntern* edge)
/**
 * @brief add an intern edge attached to the node
 * @param[in] edge pointer to the edge
 */
{
    edges_int.push_back(edge);
}

void swNode::addNeighNodes(const std::vector<swNode*> &neigh_nodes)
/**
 * @brief add neighbouring nodes
 * @param[in] neigh_nodes list of pointers to the neighbouring nodes
 */
{
    for (swNode *neigh_node : neigh_nodes)
    {
        // check if the node is already on the neighbouring nodes list
        bool found = false;
        for (swNode* node : nodes)
        {
            if (node->getSMSIdx() == neigh_node->getSMSIdx())
            {
                found = true;
            }
        }

        // if the node has not been found in the list, add it
        if (found == false)
        {
            this->nodes.push_back(neigh_node);
        }
    }
}

const std::vector<swCell*>& swNode::getCells(void)
/**
 * @brief return the list of the cell pointers
 */
{
    return cells;
}

const std::vector<swEdgeBoundary *>& swNode::getEdgesBoundary(void)
/**
 * @brief return the list of the boundary edge pointers
 */
{
    return edges_bnd;
}

const std::vector<swEdgeIntern *>& swNode::getEdgesIntern(void)
/**
 * @brief return the list of the intern edge pointers
 */
{
    return edges_int;
}

size_t swNode::getIdx(void)
/**
 * @brief return the node id (sw2d numbering style)
 */
{
    return idx_sw2d;
}

const std::vector<swNode *>& swNode::getNeighNodes(void)
/**
 * @brief return the list of the neighbouring node pointers
 */
{
    return nodes;
}

size_t swNode::getSMSIdx(void)
/**
 * @brief return the node id (SMS numbering style)
 */
{
    return idx_sms;
}

double swNode::getX(void)
/**
 * @brief return x-coordinate of the node
 */
{
    return x;
}

double swNode::getY(void)
/**
 * @brief return y-coordinate of the node
 */
{
    return y;
}

double swNode::getZ(void)
/**
 * @brief return z-coordinate of the node
 */
{
    return z;
}

void swNode::setCells(const std::vector<swCell *> cells_list)
/**
 * @brief swNode::setCells
 * @param[in] cells_list: list of cells to which current node is belonging to
 */ 
{
    cells = cells_list;
}

void swNode::setCoord(size_t sms_idx, size_t sw2d_idx, double xnod, double ynod, double znod)
/**
 * @brief swNode::setCoord
 * @param[in] sms_idx: node id in SMS (discontinuous numbering and starting from 1)
 * @param[in] sw2d_idx: node id in sw2d (continuous numbering and starting from 0)
 * @param[in] xnod: x-coordinate
 * @param[in] ynod: y-coordinate
 * @param[in] znod: z-coordinate
 */
{
    idx_sms = sms_idx;
    idx_sw2d = sw2d_idx;
    x = xnod;
    y = ynod;
    z = znod;
}