#include "triangle.h"

EquilateralTriangle::EquilateralTriangle(double sideLength, Point begin) {
    side = sideLength;
    A = begin;
    B = Point(begin.first + sideLength, begin.second);
    C = Point(begin.first + sideLength * 0.5, begin.second + sideLength * sqrt(3) * 0.5);
    center = Point(begin.first + sideLength * 0.5, begin.second + sideLength / (2 * sqrt(3)));
}

bool EquilateralTriangle::isInside(const Point& point) const {
    const double &x0 = point.first, &y0 = point.second;
    const double &x1 = A.first, &y1 = A.second;
    const double &x2 = B.first, &y2 = B.second;
    const double &x3 = C.first, &y3 = C.second;

    const double chck1 = (x1-x0)*(y2-y1)-(x2-x1)*(y1-y0);
    const double chck2 = (x2-x0)*(y3-y2)-(x3-x2)*(y2-y0);
    const double chck3 = (x3-x0)*(y1-y3)-(x1-x3)*(y3-y0);

    return !(chck1 * chck2 < 0.0 || chck1 * chck3 < 0.0 || chck2 * chck3 < 0.0);
}

double EquilateralTriangle::distanceToCenter(const Point &point) const {
    const double &x = point.first, &y = point.second;
    const double &x0 = center.first, &y0 = center.second;

    const double distance = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));

    return distance;
}

double EquilateralTriangle::getSide() const {
    return side;
}



