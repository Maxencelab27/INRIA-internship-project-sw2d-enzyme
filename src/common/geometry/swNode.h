#pragma once

#include <vector>
#include <cstddef>

class swCell;
class swEdgeBoundary;
class swEdgeIntern;

class swNode
/**
 * @class swNode
 * @brief Describe the geometrical properties of a node
 */
{
public:
    //-- General methods
    swNode(void) = default;
    swNode(size_t sms_idx, size_t sw2d_idx, double x_read, double y_read, double z_read);
    ~swNode(void) = default;

    void addCell(swCell *cell);
    void addEdgeBoundary(swEdgeBoundary *edge);
    void addEdgeIntern(swEdgeIntern *edge);
    void addNeighNodes(const std::vector<swNode*> &neigh_nodes);
    //== General methods

    //-- Getters and setters
    const std::vector<swCell *>& getCells(void);
    const std::vector<swEdgeBoundary *>& getEdgesBoundary(void);
    const std::vector<swEdgeIntern *>& getEdgesIntern(void);
    size_t getIdx(void);
    const std::vector<swNode *>& getNeighNodes(void);
    size_t getSMSIdx(void);
    double getX(void);
    double getY(void);
    double getZ(void);

    void setCells(const std::vector<swCell *> cells_list);
    void setCoord(size_t sms_idx, size_t sw2d_idx, double xnod, double ynod, double znod);
    //== Getters and setters

private:
    //-- Attributes
    std::vector<swCell*> cells;                 /**< list of the pointers to the cells attached to the node **/
    std::vector<swEdgeBoundary*> edges_bnd;     /**< list of the pointers to the boundary edges connected to the node **/
    std::vector<swEdgeIntern*> edges_int;       /**< list of the pointers to the intern edges connected to the node **/
    std::vector<swNode*> nodes;                 /**< list of the pointers to the nodes connected by an edge to the node **/
    size_t idx_sms;                              /**< node id in SMS (discontinuous numbering and starting from 1) **/
    size_t idx_sw2d;                             /**< node id in sw2d (continuous numbering and starting from 0) **/ 
    double x;                                   /**< x-coordinate of the node **/
    double y;                                   /**< y-coordinate of the node **/
    double z;                                   /**< z-coordinate of the node **/
    //== Attributes
};
