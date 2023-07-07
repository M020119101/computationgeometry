#pragma once
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

const double PI = 3.14159;

const double epsilon = 1e-10;
const double epsilonSquare = epsilon * epsilon;

#include <sstream>

template<typename T>
std::string to_string(const T& value) {
    std::ostringstream ss;
    ss << value;
    return ss.str();
}

class Matrix3D {
public:
    std::vector<std::vector<int>> a;

    Matrix3D() {
        a = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };
    }

    std::string toString() const {
        std::string result = "Matrix3D:\n";
        for (const auto& row : a) {
            result += std::to_string(row[0]) + " " + std::to_string(row[1]) + " "
                + std::to_string(row[2]) + " " + std::to_string(row[3]) + "\n";
        }
        return result;
    }

    void makeIdentical() {
        a = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };
    }

    Matrix3D multiplied(const Matrix3D& other) const {
        Matrix3D res;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                res.a[i][j] = a[i][0] * other.a[0][j]
                    + a[i][1] * other.a[1][j]
                    + a[i][2] * other.a[2][j]
                    + a[i][3] * other.a[3][j];
            }
        }
        return res;
    }

    int getDeterminant() const {
        // 计算行列式，请根据实际情况进行实现
    }

    Matrix3D getReverseMatrix() const {
        // 计算逆矩阵，请根据实际情况进行实现
    }

    //生成平移矩阵
    static Matrix3D createTranslateMatrix(int dx, int dy, int dz) {
        Matrix3D m;
        m.a[3][0] = dx;
        m.a[3][1] = dy;
        m.a[3][2] = dz;
        return m;
    }

    //缩放矩阵
    static Matrix3D createScaleMatrix(int sx, int sy, int sz) {
        Matrix3D m;
        m.a[0][0] = sx;
        m.a[1][1] = sy;
        m.a[2][2] = sz;
        return m;
    }

    //旋转矩阵
    static Matrix3D createRotateMatrix(const std::string& axis, double angle) {
        Matrix3D m;
        double sinVal = std::sin(angle);
        double cosVal = std::cos(angle);

        if (axis == "X" || axis == "x") {
            m.a[1][1] = cosVal;
            m.a[1][2] = sinVal;
            m.a[2][1] = -sinVal;
            m.a[2][2] = cosVal;
        }
        else if (axis == "Y" || axis == "y") {
            m.a[0][0] = cosVal;
            m.a[0][2] = -sinVal;
            m.a[2][0] = sinVal;
            m.a[2][2] = cosVal;
        }
        else if (axis == "Z" || axis == "z") {
            m.a[0][0] = cosVal;
            m.a[0][1] = sinVal;
            m.a[1][0] = -sinVal;
            m.a[1][1] = cosVal;
        }
        return m;
    }

    static Matrix3D createMirrorMatrix(const std::vector<int>& point, int n) {
        // 创建镜面矩阵，请根据实际情况进行实现
        Matrix3D m;
        if (n == 0) {   //0-关于XY平面镜像变换
            m.a[0][0] = -1;
            m.a[1][1] = -1;
        }
        else if (n == 1) {  //1-关于XZ平面镜像变换
            m.a[0][0] = -1;
            m.a[2][2] = -1;
        }
        else {  //2-关于YZ平面镜像变换
            m.a[1][1] = -1;
            m.a[2][2] = -1;
        }

        return m;
    }

    Matrix3D operator*(const Matrix3D& other) const {
        return multiplied(other);
    }

    Matrix3D operator+(const Matrix3D& other) const {
        Matrix3D res;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                res.a[i][j] = a[i][j] + other.a[i][j];
            }
        }
        return res;
    }

    Matrix3D operator-(const Matrix3D& other) const {
        Matrix3D res;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                res.a[i][j] = a[i][j] - other.a[i][j];
            }
        }
        return res;
    }
};


class Vector3D {
private:
    double dx, dy, dz, dw;
public:
    Vector3D(double dx = 0.0, double dy = 0.0, double dz = 0.0, double dw = 0.0) : dx(dx), dy(dy), dz(dz), dw(dw) {}

    string toString() {
        return "Vector3D: " + to_string(dx) + ", " + to_string(dy) + ", " + to_string(dz);
    }

    Vector3D clone() {
        return Vector3D(dx, dy, dz, dw);
    }

    // Getter函数
    double getDx() const { return dx; }
    double getDy() const { return dy; }
    double getDz() const { return dz; }
    double getDw() const { return dw; }

    // Setter函数
    void setDx(double val) { dx = val; }
    void setDy(double val) { dy = val; }
    void setDz(double val) { dz = val; }
    void setDw(double val) { dw = val; }

    //对当前向量取反
    void reverse() {
        dx = -dx;
        dy = -dy;
        dz = -dz;
    }

    //取反返回新向量
    Vector3D reversed() const {
        return Vector3D(-dx, -dy, -dz);
    }

    double dotProduct(Vector3D vec) {
        return dx * vec.dx + dy * vec.dy + dz * vec.dz;
    }

    Vector3D crossProduct(Vector3D vec) {
        double newX = dy * vec.dz - dz * vec.dy;
        double newY = dz * vec.dx - dx * vec.dz;
        double newZ = dx * vec.dy - dy * vec.dx;
        return Vector3D(newX, newY, newZ);
    }

    //向量放大
    void amplify(double f) {
        dx *= f;
        dy *= f;
        dz *= f;
    }

    Vector3D amplified(double f) {
        return Vector3D(dx * f, dy * f, dz * f);
    }

    double lengthSquare() {
        return dx * dx + dy * dy + dz * dz;
    }

    double length() {
        return sqrt(lengthSquare());
    }

    void normalize(){
        double len = length();
        dx /= len;
        dy /= len;
        dz /= len;
    }

    Vector3D normalized() {
        double len = length();
        return Vector3D(dx / len, dy / len, dz / len);
    }

    bool isZeroVector() {
        return lengthSquare() == 0.0;
    }

    Vector3D multiplied(Matrix3D m) {
        double x = dx * m.a[0][0] + dy * m.a[1][0] + dz * m.a[2][0] + dw * m.a[3][0];
        double y = dx * m.a[0][1] + dy * m.a[1][1] + dz * m.a[2][1] + dw * m.a[3][1];
        double z = dx * m.a[0][2] + dy * m.a[1][2] + dz * m.a[2][2] + dw * m.a[3][2];
        return Vector3D(x, y, z);
    }

    bool isParallel(Vector3D other) {
        Vector3D v = crossProduct(other);
        return v.isZeroVector();
    }

    double getAngle(Vector3D vec) {
        Vector3D v1 = normalized();
        Vector3D v2 = vec.normalized();
        double dotPro = v1.dotProduct(v2);
        if (dotPro > 1) dotPro = 1;
        else if (dotPro < -1) dotPro = -1;
        return acos(dotPro);
    }

    Vector3D getOrthoVector2D() {
        if (dx == 0) {
            return Vector3D(1, 0, 0).normalized();
        }
        else {
            return Vector3D(-dy / dx, 1, 0).normalized();
        }
    }

    //向量与x轴夹角
    double getAngle2D() {
        double rad = getAngle(Vector3D(1, 0, 0));
        if (dy < 0) {
            rad = PI * 2.0 - rad;
        }
        return rad;
    }

    Vector3D operator+(const Vector3D& other) const {
        return Vector3D(dx + other.dx, dy + other.dy, dz + other.dz);
    }

    Vector3D operator-(const Vector3D& other) const {
        return *this + other.reversed();
    }

    Vector3D operator*(double val) const {
        return Vector3D(dx * val, dy * val, dz * val);
    }
};


class Point3D {
private:
    double x, y, z, w;
public:
    Point3D(double x = 0.0, double y = 0.0, double z = 0.0, double w = 1.0) : x(x), y(y), z(z), w(w) {}

    string toString() {
        return "Point3D: " + to_string(x) + ", " + to_string(y) + ", " + to_string(z);
    }

    Point3D clone() {
        return Point3D(x, y, z, w);
    }

    // Getter函数
    double getX() const { return x; }
    double getY() const { return y; }
    double getZ() const { return z; }
    double getW() const { return w; }

    // Setter函数
    void setX(double val) { x = val; }
    void setY(double val) { y = val; }
    void setZ(double val) { z = val; }
    void setW(double val) { w = val; }

    Vector3D pointTo(Point3D other) {
        double dx = other.x - x;
        double dy = other.y - y;
        double dz = other.z - z;
        return Vector3D(dx, dy, dz);
    }

    //平移当前点
    void translate(Vector3D vec) {
        x += vec.getDx();
        y += vec.getDy();
        z += vec.getDz();
    }

    //平移后返回新点
    Point3D translated(Vector3D vec) const{
        double newX = x + vec.getDx();
        double newY = y + vec.getDy();
        double newZ = z + vec.getDz();
        return Point3D(newX, newY, newZ,vec.getDw());
    }

    Point3D multiplied(Matrix3D m) const{
        double newX = x * m.a[0][0] + y * m.a[1][0] + z * m.a[2][0] + w * m.a[3][0];
        double newY = x * m.a[0][1] + y * m.a[1][1] + z * m.a[2][1] + w * m.a[3][1];
        double newZ = x * m.a[0][2] + y * m.a[1][2] + z * m.a[2][2] + w * m.a[3][2];
        return Point3D(newX, newY, newZ);
    }

    double distance(Point3D other) {
        return pointTo(other).length();
    }

    double distanceSquare(Point3D other) {
        return pointTo(other).lengthSquare();
    }

    //两点中点
    Point3D middle(Point3D other) {
        double midX = (x + other.x) / 2;
        double midY = (y + other.y) / 2;
        double midZ = (z + other.z) / 2;
        return Point3D(midX, midY, midZ);
    }

    //重合
    bool isCoincide(Point3D other, double dis2 = epsilonSquare) {
        double d2 = pointTo(other).lengthSquare();
        if (d2 <= dis2) {
            return true;
        }
        else {
            return false;
        }
    }

    //是否完全重合
    bool isIdentical(Point3D other) {
        if (x == other.x && y == other.y && z == other.z) {
            return true;
        }
        else {
            return false;
        }
    }

    Point3D operator+(Vector3D vec) {
        return translated(vec);
    }

    Vector3D operator-(Point3D other) {
        return pointTo(other);
    }

    Point3D operator*(Matrix3D matrix) {
        double newX = x * matrix.a[0][0] + y * matrix.a[1][0] + z * matrix.a[2][0] + w * matrix.a[3][0];
        double newY = x * matrix.a[0][1] + y * matrix.a[1][1] + z * matrix.a[2][1] + w * matrix.a[3][1];
        double newZ = x * matrix.a[0][2] + y * matrix.a[1][2] + z * matrix.a[2][2] + w * matrix.a[3][2];      
        double newW = x * matrix.a[0][3] + y * matrix.a[1][3] + z * matrix.a[2][3] + w * matrix.a[3][3];
        return Point3D(newX, newY, newZ, newW);
    }
};


