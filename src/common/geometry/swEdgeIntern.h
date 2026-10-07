#pragma once

#include <vector>
#include <cstddef>

class swCell;
class swMesh;
class swNode;

class swEdgeIntern
/**
 * @class swEdgeIntern
 * @brief Describes the geometrical properties of a computational internal edge
 */
{
public:
    swEdgeIntern() = default;
    swEdgeIntern(size_t id_sw2d,
                 const std::vector<swNode *> edge_nodes,
                 const std::vector<swCell *> edge_cells);
    ~swEdgeIntern(void) = default;

    bool computeDNorm(void);
    bool correctOrientation(void);
    //== General methods

    //-- Getters and setters
    const std::vector<swCell*>& getCells();
    double getDnorm(void);
    size_t getIdx(void);
    double getLength(void);
    const std::vector<swNode*>& getNodes(void);
    //std::string getNodesId(void);
    double getNormX(void);
    double getNormY(void);
    //size_t getNS(void);
    double getXM(void);
    double getYM(void);
    double getZM(void);

    void setCells(std::vector<swCell*> &cells_list);
    void setProperties(size_t edge_idx, swNode* node_1, swNode* node_2,
                       double x, double y, double z,
                       double d, double l, double nx, double ny);

    /*void setNS(size_t id);
    void setProperties(size_t m_id,
                       const std::vector<swNode *> m_nodes,
                       double length, double Dnorm,
                       double xint, double yint, double zint);
    virtual void setNS(int id);
*/
    const swCell* getCelL();
    const swCell* getCelR();
    //== Getters and setters

private:
    //-- General methods
    bool computeCoord(void);
    bool computeLength(void);
    bool computeNormal(void);
    //== General methods

protected:
    //-- Attributes
    // geometry attributes
    double d_norm;                                  /**< distance cell-cell (or boundary-cell) in the normal direction to the edge **/
    double length;                                  /**< edge length **/
    double norm_x;                                  /**< x-component of the normal unit vector (faster than xtarray<double> for normal vector) **/
    double norm_y;                                  /**< y-component of the normal unit vector **/
    double x_m;                                       /**< x-coordinate of the edge middle point **/
    double y_m;                                       /**< y-coordinate of the edge middle point **/
    double z_m;                                       /**< z-coordinate of the edge middle point **/
    
    // other attributes
    std::vector<swCell*> cells;                     /**< list of the pointers to the cells **/
    size_t idx;                                      /**< edge index in sw2d numbering style (continuous and starting from 0) **/
    std::vector<swNode*> nodes;                     /**< list of the pointers to the edge nodes **/
    //== Attributes
};
