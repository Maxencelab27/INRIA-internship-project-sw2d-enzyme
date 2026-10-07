#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

#include "swGeometry.h"

#include "../common/swLog.h"

double swGeometry::computeArea(std::vector<double> x_poly, std::vector<double> y_poly)
/**
 * @brief Compute the area of a polygon from the coordinates of its vertices
 * @param[in] x_poly x coordinates of the vertices
 * @param[in] y_poly y coordinates of the vertices
 * @return area of the polygon 
 */
{
    // Sort vertices clockwise
    sortVertices(x_poly, y_poly);

    // Compute the area : sum of the cross product of each pair of consecutive vertices
    double area = 0;
    size_t n_vert = x_poly.size();
    for (size_t k=0; k<n_vert; k++)
    {
        double cp = crossProduct(x_poly[k], y_poly[k], x_poly[(k+1)%n_vert], y_poly[(k+1)%n_vert]);
        area += cp;
    }
    area /= 2;
    return (area);
}

double swGeometry::computeIntersectionArea(std::vector<double> x_poly_1, std::vector<double> y_poly_1,
                                           std::vector<double> x_poly_2, std::vector<double> y_poly_2)
/**
 * @brief compute the area of the intersection between two polygons
 * @param[in] x_poly_1 x coordinates of the first polygone
 * @param[in] y_poly_1 y coordinates of the first polygon
 * @param[in] x_poly_2 x coordinates of the second polygon
 * @param[in] y_poly_2 y coordinates of the second polygon
 * @return double area of intersection 
 */
{
    // Sort the vertices of poly_1 and poly_2;
    sortVertices(x_poly_1, y_poly_1);
    sortVertices(x_poly_2, y_poly_2);
    
    //-- Compute and store the coordinates of the intersection polygon (poly_3)
    std::vector<double> x_poly_3;
    std::vector<double> y_poly_3;

    size_t n_vert_1 = x_poly_1.size();
    size_t n_vert_2 = x_poly_2.size();

    // Add every vertex of poly_1 contained in poly_2
    for (size_t k = 0; k < n_vert_1; k++) 
    {
        if ((isInsidePolygon(x_poly_1[k], y_poly_1[k], x_poly_2, y_poly_2, 0)))
        {
            x_poly_3.push_back(x_poly_1[k]);
            y_poly_3.push_back(y_poly_1[k]);
        }
    }

    // Add every vertex of poly_2 contained in poly_1
    for (size_t k = 0; k < n_vert_2; k++) 
    {
        if ((isInsidePolygon(x_poly_2[k], y_poly_2[k], x_poly_1, y_poly_1, 0)))
        {
            x_poly_3.push_back(x_poly_2[k]);
            y_poly_3.push_back(y_poly_2[k]);
        }
    }

    // Add every intersection points between sides of poly_1 and sides of poly_2
    for (size_t k_1 = 0; k_1 < n_vert_1; k_1++)
    {
        std::vector<double> x_side_1 = {x_poly_1[k_1], x_poly_1[(k_1+1) % n_vert_1]};
        std::vector<double> y_side_1 = {y_poly_1[k_1], y_poly_1[(k_1+1) % n_vert_1]};
        for (size_t k_2 = 0; k_2 < n_vert_2; k_2++)
        {
            std::vector<double> x_side_2 = {x_poly_2[k_2], x_poly_2[(k_2+1) % n_vert_2]};
            std::vector<double> y_side_2 = {y_poly_2[k_2], y_poly_2[(k_2+1) % n_vert_2]};
            std::vector<double> intersection_pt = {0,0};

            if ((segmentsIntersection(x_side_1, y_side_1, x_side_2, y_side_2, intersection_pt)) && (!isVertex(x_poly_3, y_poly_3, intersection_pt[0], intersection_pt[1])))
            {
                x_poly_3.push_back(intersection_pt[0]);
                y_poly_3.push_back(intersection_pt[1]);
            }
        }
    }
    //== Compute and store the coordinates of the intersection polygon (poly_3)

    // Compute and return area of poly_3
    return (computeArea(x_poly_3, y_poly_3));
}

double swGeometry::crossProduct(double x1, double y1, double x2, double y2)
/**
 * @brief compute the cross product bewteen two vectors
 * @param[in] x1 x-coordinate of the first vector
 * @param[in] y1 y-coordinate of the first vector
 * @param[in] x2 x-coordinate of the second vector
 * @param[in] y2 t-coordinate of the second vector
 * @return double cross product result
 */
{
    return (x1*y2 - x2*y1);
}

bool swGeometry::isInsidePolygon(double const& x, double const& y, std::vector<double> x_poly, std::vector<double> y_poly, double small)
/**
 * @brief check if a point is inside a polygon
 * @param[in] x x coordinate of the point
 * @param[in] y y coordinate of the point
 * @param[in] x_poly x coordinates of the polygon
 * @param[in] y_poly y coordinates of the polygon
 * @param[in] small tolerance to test negative values
 * @return boolean indicating if the point is inside the polygon or not
 */
{
    sortVertices(x_poly, y_poly);
    
    int intersections = 0;

    size_t n_sides = x_poly.size();
    for (size_t k = 0; k < n_sides; k++) {
        double x1 = x_poly[k];
        double y1 = y_poly[k];
        double x2 = x_poly[(k + 1) % n_sides];
        double y2 = y_poly[(k + 1) % n_sides];

        if (std::abs(y1-y2) < small) {continue ;} // Side is horizontal -> paralelle with horizontal line

        double x_intersection = (x2 - x1) * (y - y1) / (y2 - y1) + x1;
        if ( (x_intersection > x) && ( ((y1 < y) && y2 >= y) || (y2 < y && y1 >= y)) )
        {
            intersections++;
        }
    }
    return (intersections % 2 == 1);
}

bool swGeometry::isOnPolygonContour(double const& x, double const& y, std::vector<double> x_poly, std::vector<double> y_poly)
/**
 * @brief check if a point is on the contour of a polygon
 * @param[in] x x coordinate of the point
 * @param[in] y y coordinate of the point
 * @param[in] x_poly x coordinates of the polygon
 * @param[in] y_poly y coordinates of the polygon
 * @return boolean indicating if the point is on the polygon contour or not
 */
{
    swGeometry::sortVertices(x_poly, y_poly);
    
    size_t n_sides = x_poly.size();
    for (size_t k=0; k<n_sides; k++)
    {
        std::vector<double> x_side = {x_poly[k],x_poly[(k+1) % n_sides]};
        std::vector<double> y_side = {y_poly[k],y_poly[(k+1) % n_sides]};

        if (isOnSegment(x, y, x_side, y_side))
        {
            return(true);
        }
    }
    return(false);
}

bool swGeometry::isOnSegment(double const& x, double const& y, std::vector<double> const& x_segment, std::vector<double> const& y_segment)
/**
 * @brief check if a point is on a segment
 * @param[in] x x coordinate of the point
 * @param[in] y y coordinate of the point
 * @param[in] x_segment x coordinates of the segment
 * @param[in] y_segment y coordinates of the segment
 * @return boolean indicating if the point is on the segment or not
 */
{    
    if (y_segment[1] == y_segment[0])
    {
        return((y == y_segment[0]) &&
               (x >= std::min(x_segment[0], x_segment[1])) &&
               (x <= std::max(x_segment[0], x_segment[1])));
    }

    else
    {
        return ((x == (y-y_segment[0]) * (x_segment[1]-x_segment[0]) / (y_segment[1]-y_segment[0]) + x_segment[0]) &&
                (y >= std::min(y_segment[0], y_segment[1])) &&
                (y <= std::max(y_segment[0], y_segment[1])));
    }
}

bool swGeometry::segmentsIntersection(std::vector<double> const& x_segment_1, std::vector<double> const& y_segment_1,
                                        std::vector<double> const& x_segment_2, std::vector<double> const& y_segment_2, std::vector<double>& point)
/**
 * @brief Check if two segment intersect and compute the intersection point if so
 * @param[in] x_segment_1 x coordinates of the first segment
 * @param[in] y_segment_1 y coordinates of the first segment
 * @param[in] x_segment_2 x coordinates of the second segment
 * @param[in] y_segment_2 y coordinates of the second segment
 * @param[in,out] point point to which the value of the intersection point is assigned
 * @return boolean indicating if the segments intersect or not
 */
{
    //-- Compute the intersection point of both lines
    double x_inter, y_inter;
    
    if ((x_segment_1[0] != x_segment_1[1]) && (x_segment_2[0] != x_segment_2[1])) // none of the segments are vertical
    {
        double a_1 = (y_segment_1[1] - y_segment_1[0]) / (x_segment_1[1] - x_segment_1[0]);
        double b_1 = y_segment_1[0] - a_1*x_segment_1[0];

        double a_2 = (y_segment_2[1] - y_segment_2[0]) / (x_segment_2[1] - x_segment_2[0]);
        double b_2 = y_segment_2[0] - a_2*x_segment_2[0];

        if (a_1 == a_2) // lines are paralleles
        {
            return (false);
        }

        x_inter = - (b_2 - b_1) / (a_2 - a_1);
        y_inter = a_1 * x_inter + b_1; 
    }
    else if ((x_segment_1[0] == x_segment_1[1]) && (x_segment_2[0] != x_segment_2[1])) // only segment 1 is vertical
    {
        double a_2 = (y_segment_2[1] - y_segment_2[0]) / (x_segment_2[1] - x_segment_2[0]);
        double b_2 = y_segment_2[0] - a_2*x_segment_2[0];
        
        x_inter = x_segment_1[0];
        y_inter = a_2 * x_inter + b_2;
    }
    else if ((x_segment_1[0] != x_segment_1[1]) && (x_segment_2[0] == x_segment_2[1])) // only segment 2 is vertical
    {
        double a_1 = (y_segment_1[1] - y_segment_1[0]) / (x_segment_1[1] - x_segment_1[0]);
        double b_1 = y_segment_1[0] - a_1*x_segment_1[0];

        x_inter = x_segment_2[0];
        y_inter = a_1 * x_inter + b_1;
    }
    else // both segments are vertical
    {
        return(false);
    }
    
    //== Compute the intersection point of both lines

    //-- Check if the point is inside both segments
    if (((std::min(x_segment_1[0], x_segment_1[1]) <= x_inter) && (x_inter <= std::max(x_segment_1[0], x_segment_1[1]))) &&
        ((std::min(x_segment_2[0], x_segment_2[1]) <= x_inter) && (x_inter <= std::max(x_segment_2[0], x_segment_2[1]))) &&
        ((std::min(y_segment_1[0], y_segment_1[1]) <= y_inter) && (y_inter <= std::max(y_segment_1[0], y_segment_1[1]))) &&
        ((std::min(y_segment_2[0], y_segment_2[1]) <= y_inter) && (y_inter <= std::max(y_segment_2[0], y_segment_2[1]))))
    {
        point[0] = x_inter;
        point[1] = y_inter;
        return(true);
    }
    else
    {
        return(false);
    }
    //== Check if the point is inside both segments
}

bool swGeometry::isVertex(std::vector<double> const& x_poly, std::vector<double> const& y_poly, double const& x, double const& y)
/**
 * @brief Check if a point belongs to a polygon vertices
 * @param[in] x_poly x coordinates of the polygon
 * @param[in] y_poly y coordinates of the polygon
 * @param[in] x x coordinate of the point
 * @param[in] y y coordinate of the point
 * @return boolean indicating if the point is a vertex of the polygon or not
 */
{
    for(size_t k = 0; k < x_poly.size(); k++)
    {
        if (x_poly[k]==x && y_poly[k]==y)
        {
            return (true);
        }
    }
    return(false);
}

double swGeometry::lenghtSegmentInPolygon(std::vector<double> x_poly, std::vector<double> y_poly, std::vector<double> const& x_segment, std::vector<double> const& y_segment)
/**
 * @brief Compute the lenght of intersection between a segment and a polygon
 * @param[in] x_poly x coordinates of the polygon
 * @param[in] y_poly y coordinates of the polygon
 * @param[in] x_segment x coordinates of the segment
 * @param[in] y_segment y coordinates of the segment
 * @return lenght of intersection
 */
{
    ///-- Get the coordinates of every point of intersection between the segment and the side of the polygon
    std::vector<double> x_intersections {x_segment[0], x_segment[1]};
    std::vector<double> y_intersections {y_segment[0], y_segment[1]};

    swGeometry::sortVertices(x_poly, y_poly);

    size_t n_sides = x_poly.size();
    for (size_t k=0; k<n_sides; k++)
    {
        std::vector<double> x_side = {x_poly[k], x_poly[(k+1) % n_sides]};
        std::vector<double> y_side = {y_poly[k], y_poly[(k+1) % n_sides]};

        std::vector<double> point_inter = {0,0};
        if (segmentsIntersection(x_side, y_side, x_segment, y_segment, point_inter))  // check if the segment intersects the side
        {
            if(!isVertex(x_intersections, y_intersections, point_inter[0], point_inter[1])) // check that the intersection point is not already in the list (vertices)
            {
                x_intersections.push_back(point_inter[0]);
                y_intersections.push_back(point_inter[1]);
            }
        }
    }
    //== Get the coordinates of every point of intersection between the segment and the side of the polygon

    ///-- Sort intersection points by x increasing (y if x constant)
    size_t n_points = x_intersections.size();

    if (x_intersections[0] != x_intersections[1])
    {
        for (size_t i=0; i<n_points; i++)
        {
            size_t j_min = i;
            for (size_t j=i; j<n_points; j++)
            {
                if (x_intersections[j] < x_intersections[j_min])
                {
                    j_min = j;
                }
            }
            std::swap(x_intersections[i],x_intersections[j_min]);
            std::swap(y_intersections[i],y_intersections[j_min]);
        }
    }
    else // x constant -> sort by y increasing
    {
        for (size_t i=0; i<n_points; i++)
        {
            size_t j_min = i;
            for (size_t j=i; j<n_points; j++)
            {
                if (y_intersections[j] < y_intersections[j_min])
                {
                    j_min = j;
                }
            }
            std::swap(x_intersections[i],x_intersections[j_min]);
            std::swap(y_intersections[i],y_intersections[j_min]);
        }
    }
    //== Sort intersection points by x increasing (y if x constant)

    ///-- Check if the first portion of the segment is inside the polygon
    bool inside;

    if (isOnPolygonContour(x_intersections[0], y_intersections[0], x_poly, y_poly)) // if the first point is on an side, check the middle point of the first portion
    {
        double x_middle = (x_intersections[0] + x_intersections[1]) / 2;
        double y_middle = (y_intersections[0] + y_intersections[1]) / 2;
        inside = swGeometry::isInsidePolygon(x_middle, y_middle, x_poly, y_poly, 0) || isOnPolygonContour(x_middle, y_middle, x_poly, y_poly);
    }
    else
    {
        inside = swGeometry::isInsidePolygon(x_intersections[0], y_intersections[0], x_poly, y_poly, 0);
    }
    //== Check if the first portion of the segment is inside the polygon

    ///-- Compute and add the lenght of each portion of the segment inside the polygone
    double lenght = 0;
    for (size_t k=0; k<x_intersections.size()-1; k++)
    {
        if (inside)
        {
            double x = x_intersections[k+1] - x_intersections[k];
            double y = y_intersections[k+1] - y_intersections[k];
            lenght += std::sqrt(x*x + y*y);
        }
        inside = !inside;
    }
    //== Compute and add the lenght of each portion of the segment inside the polygone

    return(lenght);
}

std::vector<double> swGeometry::rotatePoint(double x, double y, double angle, std::vector<double> point)
/**
 * @brief Compute new coordinates of a given point after a rotation by a given angle degree around a given point
 * @param[in] x x coordinate of original point
 * @param[in] y y coordinate of original point
 * @param[in] angle rotation angle (in degrees)
 * @param[in] point rotation point
 * @return std::vector<double> the coordinates of the rotated point
 */
{
    // convert degree angle into radians
    double angle_rad = angle*M_PI/180;

    // center of rotation coordinates
    double x_c = point[0];
    double y_c = point[1];

    // new rotated point coordinates
    double x_rot = (x-x_c)*cos(angle_rad) - (y-y_c)*sin(angle_rad) + x_c;
    double y_rot = (x-x_c)*sin(angle_rad) + (y-y_c)*cos(angle_rad) + y_c;

    std::vector<double> rot_coord;
    rot_coord.push_back(x_rot);
    rot_coord.push_back(y_rot);

    return (rot_coord);
}

bool swGeometry::sortVertices(std::vector<double>& x_poly, std::vector<double>& y_poly)
/**
 * @brief Sort the vertices of a polygon clockwise around the center of the polygon
 * @param[in,out] x_poly x coordinates of the polygon
 * @param[in,out] y_poly coordinates of the polygon
 * @return EXIT_SUCCESS/EXIT_FAILURE 
 */
{    
    /// Compute the center of the polygon (average coordinates)
    double x_c = 0;
    double y_c = 0;
    for (size_t k=0; k<x_poly.size(); k++)
    {
        x_c += x_poly[k];
        y_c += y_poly[k];
    }
    x_c /= x_poly.size();
    y_c /= y_poly.size();

    /// Compute the angles between(x-i-x_c, y_i-y_c) and (1, 0) vectors
    std::vector<double> angles;
    for (size_t k=0; k<x_poly.size(); k++)
    {
        double angle = atan2(y_poly[k]-y_c,x_poly[k]-x_c);
        if (angle < 0)
        {
            angle += 2*M_PI;
        }
        angles.push_back(angle);
    }

    /// Sort points by angle
    size_t lenght = x_poly.size();
    for (size_t i=0; i<lenght; i++)
    {
        size_t j_min = i;
        for (size_t j=i; j<lenght; j++)
        {
            if (angles[j] < angles[j_min])
            {
                j_min = j;
            }
        }

        std::swap(x_poly[i],x_poly[j_min]);
        std::swap(y_poly[i],y_poly[j_min]);
        std::swap(angles[i],angles[j_min]);
    }

    return(EXIT_SUCCESS);
}

bool swGeometry::triangleOrientation(std::vector<double> const& x_poly, std::vector<double> const& y_poly, size_t const& a, size_t const& b, size_t const& c)
/**
 * @brief check the orientation of the triangle formed by three vertices of a polygon
 * @param[in] x_poly x coordinates of the polygon
 * @param[in] y_poly y coordinates of the polygon
 * @param[in] a index of the first point in the coordinates lists
 * @param[in] b index of the second point in the coordinates lists
 * @param[in] c index of the third point in the coordinates lists
 * @return bool: true if the points are clockwise, else false
 */
{
    return ((crossProduct(x_poly[a], y_poly[a], x_poly[b], y_poly[b]) +
             crossProduct(x_poly[b], y_poly[b], x_poly[c], y_poly[c]) +
             crossProduct(x_poly[c], y_poly[c], x_poly[a], y_poly[a])) < 0);
}
