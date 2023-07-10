#pragma once
#include <iostream>
#include <cmath>
#include "GeomBase.h"

class Ray {
public:
    Point3D P;
    Vector3D V;

    Ray(const Point3D& P, Vector3D& V) : P(P), V(V) {
        V.normalize();
    }

    std::string toString() const {
        return "Ray\nP " + pointToString(P) + "\nV " + vectorToString(V) + "\n";
    }

private:
    std::string pointToString(const Point3D& point) const {
        return "(" + std::to_string(point.getX()) + ", " + std::to_string(point.getY()) + ", " + std::to_string(point.getZ()) + ")";
    }

    std::string vectorToString(const Vector3D& vector) const {
        return "<" + std::to_string(vector.getDx()) + ", " + std::to_string(vector.getDy()) + ", " + std::to_string(vector.getDz()) + ">";
    }
};
