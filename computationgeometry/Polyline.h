#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include "GeomBase.h"
#include <tuple>
#include "Segment.h"

class Polyline {
private:
    std::vector<Point3D> points;

public:
    Polyline() {}

    std::size_t count() const {
        return points.size();
    }

    void addPoint(const Point3D& pt) {
        points.push_back(pt);
    }

    const Point3D& startPoint() const {
        return points.front();
    }

    const Point3D& endPoint() const {
        return points.back();
    }

    const std::vector<Point3D> getpoints() const {
        return points;
    }

    void addTuple(const std::tuple<double, double, double>& tuple) {
        double x = std::get<0>(tuple);
        double y = std::get<1>(tuple);
        double z = std::get<2>(tuple);
        points.push_back(Point3D(x, y, z));
    }

    //开头加入一个点
    void raddPoint(const Point3D& pt) {
        points.insert(points.begin(), pt);
    }

    //根据点的序号删除一个点
    std::vector<Point3D> removePoint(int index) {
        //Point3D removedPoint = points[index];
        points.erase(points.begin() + index);
        return points;
    }

    //根据点的序号获取一个点
    Point3D point(int index) const {
        return points[index];
    }

    bool isClosed(){
        if (count() <= 2) {
            return false;
        }
        return points.front().isCoincide(points.back());
    }

    double getArea() const {
        double area = 0.0;
        for (int i = 0; i < count() - 1; i++) {
            area += 0.5 * (points[i].getX() * points[i + 1].getY() - points[i + 1].getX() * points[i].getY());
        }
        return area;
    }

    void reverse() {
        int sz = count();
        for (int i = 0; i < sz / 2; i++) {
            std::swap(points[i], points[sz - 1 - i]);
        }
    }

    //逆时针
    void makeCCW() {
        if (getArea() < 0) {
            reverse();
        }
    }

    void makeCW() {
        if (getArea() > 0) {
            reverse();
        }
    }

    bool isCCW() {
        return getArea() > 0;
    }

    //平移
    void translate(const Vector3D& vec) {
        for (int i = 0; i < points.size(); i++) {
            points[i].translate(vec);
        }
    }

    //添加线段
    bool appendSegment(Segment& seg) {
        if (count() == 0) {
            points.push_back(seg.A);
            points.push_back(seg.B);
        }
        else {
            Point3D startPnt = startPoint();
            Point3D endPnt = endPoint();

            if (seg.A.isCoincide(endPnt)) {
                addPoint(seg.B);
            }
            else if (seg.B.isCoincide(endPnt)) {
                addPoint(seg.A);
            }
            else if (seg.A.isCoincide(startPnt)) {
                raddPoint(seg.B);
            }
            else if (seg.B.isCoincide(startPnt)) {
                raddPoint(seg.A);
            }
            else {
                return false;
            }
        }

        return true;
    }
};

void writePolyline(const std::string& path, const Polyline& polyline) {
    std::ofstream file(path);
    if (file.is_open()) {
        file << polyline.count() << '\n';
        for (const Point3D& pt : polyline.getpoints()) {
            file << pt.getX() << ',' << pt.getY() << ',' << pt.getZ() << '\n';
        }
        file.close();
    }
    else {
        std::cout << "Failed to open file." << std::endl;
    }
}

Polyline readPolyline(const std::string& path) {
    std::ifstream file(path);
    Polyline poly;

    if (file.is_open()) {
        std::size_t number;
        file >> number;

        for (std::size_t i = 0; i < number; ++i) {
            double x, y, z;
            char comma;
            file >> x >> comma >> y >> comma >> z;
            poly.addPoint(Point3D(x, y, z));
        }

        file.close();
    }
    else {
        std::cout << "Failed to open file." << std::endl;
    }

    return poly;
}
