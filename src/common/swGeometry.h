#pragma once

#include <vector>

namespace swGeometry
/**
* @brief namespace containing geometry functions
*/
{
    ///-- General methods
    double computeArea(std::vector<double> x_poly, std::vector<double> y_poly);
    double computeIntersectionArea(std::vector<double> x_poly_1, std::vector<double> y_poly_1,
                                    std::vector<double> x_poly_2, std::vector<double> y_poly_2);
    double crossProduct(double x1, double y1, double x2, double y2);
    bool isInsidePolygon(double const& x, double const& y, std::vector<double> x_poly, std::vector<double> y_poly, double small);
    bool isOnPolygonContour(double const& x, double const& y, std::vector<double> x_poly, std::vector<double> y_poly);
    bool isOnSegment(double const& x, double const& y, std::vector<double> const& x_segment, std::vector<double> const& y_segment);
    bool isVertex(std::vector<double> const& x_poly, std::vector<double> const& y_poly, double const& x, double const& y);
    double lenghtSegmentInPolygon(std::vector<double> x_poly, std::vector<double> y_poly, std::vector<double> const& x_segment, std::vector<double> const& y_segment);
    std::vector<double> rotatePoint(double x_coord, double y_coord, double angle, std::vector<double> point);
    bool segmentsIntersection(std::vector<double> const& x_segment_1, std::vector<double> const& y_segment_1,
                                std::vector<double> const& x_segment_2, std::vector<double> const& y_segment_2, std::vector<double>& point);
    bool sortVertices(std::vector<double>& x_poly, std::vector<double>& y_poly);
    bool triangleOrientation(std::vector<double> const& x_poly, std::vector<double> const& y_poly, size_t const& a, size_t const& b, size_t const& c);
    //== General methods
}
