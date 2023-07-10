#pragma once
#include <iostream>
#include <cmath>
#include "GeomBase.h"

class Segment {
public:
    Point3D A, B;

    Segment(const Point3D& A, const Point3D& B) : A(A), B(B) {}

    std::string toString() const {
        return "Segment\nA " + pointToString(A) + "\nB " + pointToString(B) + "\n";
    }

    double length() const {
        return distance(A, B);
    }

    Point3D direction() const {
        return pointTo(B, A);
    }

    void swap() {
        std::swap(A, B);
    }

    void multiply(double m) {
        A = multiplied(A, m);
        B = multiplied(B, m);
    }

    Segment multiplied(double m) const {
        return Segment(multiplied(A, m), multiplied(B, m));
    }

private:
    std::string pointToString(const Point3D& point) const {
        return "(" + std::to_string(point.getX()) + ", " + std::to_string(point.getY()) + ", " + std::to_string(point.getZ()) + ")";
    }

    double distance(const Point3D& p1, const Point3D& p2) const {
        double dx = p2.getX() - p1.getX();
        double dy = p2.getY() - p1.getY();
        double dz = p2.getZ() - p1.getZ();
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    Point3D pointTo(const Point3D& p1, const Point3D& p2) const {
        return Point3D(p2.getX() - p1.getX(), p2.getY() - p1.getY(), p2.getZ() - p1.getZ());
    }

    Point3D multiplied(const Point3D& point, double m) const {
        return Point3D(point.getX() * m, point.getY() * m, point.getZ() * m);
    }
};

