#pragma once

#include <algorithm>
#include <vector>
#include <string>

class swEdgeBoundary;
class swEdgeIntern;
class swNode;


class swCell
/**
 * @class swCell
 * @brief Describe the geometrical properties of a computational cell
 */
{
public:
    //-- General methods
    swCell() = default;
    swCell(size_t i_id_sms, size_t i_id_sw2d, std::vector<swNode*> nodes_ptr);
    ~swCell(void) = default;

    void addEdgeBoundary(swEdgeBoundary *edge);
    void addEdgeIntern(swEdgeIntern *edge);
    void addNeigh(swCell *cell);
    bool computeArea(void);
    bool computeGravityCenter(void);
    bool computePerimeter(void);
    bool correctTopography(std::vector<double> corr_zb, const std::string &param_corr);
    //== General methods

    //-- Getters and setters
    double getArea(void);
    const std::vector<swEdgeBoundary*>& getEdgesBoundary(void);
    const std::vector<swEdgeIntern*>& getEdgesIntern(void);
    size_t getIdx(void);
    const std::vector<swCell*>& getNeigh(void);
    const std::vector<swNode*>& getNodes(void);
    double getPerimeter(void);
    size_t getSMSIdx(void);
    double getXG(void);
    double getYG(void);
    double getZBottom(void);

    void setEdgesBnd(const std::vector<swEdgeBoundary *> &edges_bnd_list);
    void setEdgesInt(const std::vector<swEdgeIntern *> &edges_int_list);
    void setNeigh(const std::vector<swCell *> &neigh_cells_list);
    void setNodes(const std::vector<swNode *> &nodes_list);
    void setProperties(size_t idx_sw2d, size_t sms_id,
                       double x, double y, double z,
                       double a, double p,
                       size_t nb_n, size_t nb_bnd, size_t nb_int, size_t nb_nei);
    //== Getters and setters

private:
    //-- Attributes
    // General attributes
    size_t id_sms = 0;                      /**< cell id in SMS (discontinuous numbering and starting from 1)**/
    size_t idx = 0;                     /**< cell index in sw2d (continuous numbering and starting from 0) **/

    // Elements linked to the cell
    std::vector<swEdgeBoundary*> edges_bnd; /**< list of the pointers to the cell boun edges **/
    std::vector<swEdgeIntern*> edges_int;   /**< list of the pointers to the cell intern edges **/
    std::vector<swCell*> neigh_cells;       /**< list of the pointers to the neighbouring cells **/
    std::vector<swNode*> nodes;             /**< list of the pointers to the cell nodes **/

    // Cell geometry
    double area = 0;                        /**< cell area **/
    double perimeter = 0;                   /**< cell perimeter **/
    double x_g = 0;                          /**< x-coordinate of the cell gravity center **/
    double y_g = 0;                          /**< y-coordinate of the cell gravity center **/
    double z_b = 0;                          /**< bottom elevation of the cell **/
    //== Attributes
};
