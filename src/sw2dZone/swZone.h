#pragma once

#include "../common/geometry/swCell.h"
#include "../common/geometry/swEdgeBoundary.h"
#include "../common/geometry/swEdgeIntern.h"
#include "../common/geometry/swMesh.h"

class swZone
{
public:
    //-- General methods
    swZone(swMesh *mesh, std::string name, std::vector<double>& m_x_list, std::vector<double>& m_y_list);
    ~swZone() = default;

    bool addVariable(std::string var, double value, std::vector<std::string> default_var);
    bool assignElements();
    static swZone *create(swMesh *mesh, std::string &line);
    bool isInside(double x, double y);
    //== General methods

    //-- Getters
    std::vector<swCell*> getCellList();
    std::vector<swEdgeBoundary *> getEdgeBoundaryList();
    std::vector<swEdgeIntern *> getEdgeInternList();
    std::string getName();
    double getValue(size_t var_id);
    std::vector<std::string> getVariables();
    //== Getters

private:
    //-- Attributes
    std::vector<swCell*> cell_list;
    std::vector<swEdgeBoundary *> edge_bnd_list;
    std::vector<swEdgeIntern *> edge_int_list;
    swMesh *m_mesh = nullptr;
    std::string m_name;
    std::vector<double> values;
    std::vector<std::string> variables;
    std::vector<double> x_list;
    std::vector<double> y_list;
    //== Attributes
};


