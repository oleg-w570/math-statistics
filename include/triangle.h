#pragma once
#include <utility>
#include <cmath>

using Point = std::pair<double, double>;

class EquilateralTriangle {
    double side;
    Point A;
    Point B;
    Point C;
    Point center;

public:
    explicit EquilateralTriangle(double sideLength = 1.0, Point begin = Point(0.0, 0.0));
    bool isInside(const Point& point) const;
    double distanceToCenter(const Point& point) const;
    double getSide() const;
};