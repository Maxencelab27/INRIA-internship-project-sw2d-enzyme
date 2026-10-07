#pragma once

#include <vector>
#include <cstddef>

class swBoundaryCondition;
class swCell;
class swNode;

class swEdgeBoundary
/**
 * @class swEdge
 * @brief Describes the geometrical properties of a computational boundary edge
 */
{
public:
    //-- General methods
    swEdgeBoundary() = default;
    swEdgeBoundary(size_t id_sw2d,
                   const std::vector<swNode *> edge_nodes,
                   swCell *edge_cell,
                   int code);
    ~swEdgeBoundary(void) = default;

    bool computeDNorm(void);  
    bool correctOrientation(void);
    //== General methods

    //-- Getters and setters
    swCell *getCell();
    double getDnorm(void);
    size_t getIdx(void);
    double getLength(void);
    const std::vector<swNode*>& getNodes(void);
    double getNormX(void);
    double getNormY(void);
    double getXM(void);
    double getYM(void);
    double getZM(void);

    void setCell(swCell *cell_pointer);
    void setProperties(size_t edge_idx, swNode* node_1, swNode* node_2,
                       double x, double y, double z,
                       double d, double l, double nx, double ny,
                       int code);

    int getBCCode(void);
    void setBCCode(int code);
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
    double d_norm;                                  /**< distance boundary-cell in the normal direction to the edge **/
    double length;                                  /**< edge length **/
    double norm_x;                                  /**< x-component of the normal unit vector (faster than xtarray<double> for normal vector) **/
    double norm_y;                                  /**< y-component of the normal unit vector **/
    double x_m;                                       /**< x-coordinate of the edge middle point **/
    double y_m;                                       /**< y-coordinate of the edge middle point **/
    double z_m;                                       /**< z-coordinate of the edge middle point **/
    
    // other attributes
    swCell *cell;                                   /**< pointers to the cell **/
    size_t idx;                                     /**< edge index in sw2d numbering style (continuous and starting from 0) **/
    std::vector<swNode*> nodes;                     /**< list of the pointers to the edge nodes **/

    swBoundaryCondition *m_bc = nullptr;
    int bc_code = 0;                               /**< boundary condition code */
    //== Attributes
};
